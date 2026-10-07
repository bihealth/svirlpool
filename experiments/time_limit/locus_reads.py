"""Step 1 of the read-selection study: the alignment structure of slow loci.

For each studied container and each of its CRs, the reads the consensus stage
fetches (svirlpool's ReadSequenceCache: primary + supplementary alignments
overlapping the CR) and, per read, how it crosses the CR:

  single    one alignment covers the whole CR
  split     the read crosses the CR in >= 2 alignments of the same strand,
            collinear on the read (an allele longer than, or diverged from,
            the reference: e.g. an expanded satellite)
  left / right / inside   anchored on one side only, or not at all
  anchor    for single / split: min(CR start - start of the left alignment,
            end of the right alignment - CR end), i.e. how far the read runs
            continuously beyond the CR on its shorter side

Per container (locus features): alignments, reads, supplementary share,
median aligned length of the alignments at the CR, reads crossing it (single
/ split), their anchors, and the share of the alignments that are short
fragments (< 2 kbp aligned) of reads that do not cross the CR.

Writes <out>/containers.tsv, <out>/reads.tsv and <out>/reads.fa (the full
sequences of every crossing read and of up to --sample other reads per CR),
for the truth step (locus_truth.py).

usage: locus_reads.py <svp root> <variant> <out dir> [--min-seconds 120]
         [--controls 100] [--extra <crs_containers.db>:<crID>,...]
"""

from __future__ import annotations

import argparse
import random
import re
import sys
from pathlib import Path

import container_value as CV
import numpy as np
import pandas as pd
import pysam

from svirlpool.localassembly import consensus as C
from svirlpool.localassembly import read_cache as RC

BAM = "/home/mayv_c/biodata/local/alignments/HG002.minimap2.20x.softclipped.bam"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
NEAR = 200_000  # alignments of a read this far from the CR are not its anchors


def sa_alignments(aln: pysam.AlignedSegment) -> list[dict]:
    """The read's alignments (this one plus its SA tag) as dicts with ref and
    read coordinates; read coordinates on the original read orientation."""
    out = []
    qlen = aln.infer_read_length()
    entries = [(aln.reference_name, aln.reference_start + 1, "-" if aln.is_reverse else "+",
                aln.cigarstring, aln.mapping_quality)]
    if aln.has_tag("SA"):
        for e in aln.get_tag("SA").rstrip(";").split(";"):
            f = e.split(",")
            entries.append((f[0], int(f[1]), f[2], f[3], int(f[4])))
    seen = set()
    for rname, pos, strand, cigar, mapq in entries:
        if (rname, pos, strand) in seen:
            continue
        seen.add((rname, pos, strand))
        ops = [(int(n), o) for n, o in CIG.findall(cigar)]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        trail = ops[-1][0] if ops[-1][1] in "SH" else 0
        rlen = sum(n for n, o in ops if o in "MDN=X")
        qs, qe = lead, qlen - trail
        if strand == "-":
            qs, qe = qlen - qe, qlen - qs
        out.append({"chr": rname, "rs": pos - 1, "re": pos - 1 + rlen, "strand": strand,
                    "qs": qs, "qe": qe, "mapq": mapq, "ops": ops, "lead": lead, "qlen": qlen})
    return out


def ref_to_read(a: dict, x: int) -> int:
    """Read coordinate (original read orientation) of reference position x
    through alignment a (x clamped to the alignment)."""
    x = min(max(x, a["rs"]), a["re"])
    r, q = a["rs"], a["lead"]
    for n, o in a["ops"]:
        if o in "M=X":
            if r + n >= x:
                q += x - r
                break
            r += n
            q += n
        elif o in "DN":
            if r + n >= x:
                break
            r += n
        elif o == "I":
            q += n
    return a["qlen"] - q if a["strand"] == "-" else q


def read_interval(a: dict, b: dict, s: int, e: int) -> tuple[int, int]:
    """The read's segment from reference position s (on a) to e (on b)."""
    q1, q2 = ref_to_read(a, s), ref_to_read(b, e)
    return min(q1, q2), max(q1, q2)


def crossing(alns: list[dict], chrom: str, s: int, e: int) -> tuple[str, int, int, int]:
    """How a read crosses [s, e): (kind, anchor, q1, q2), q1-q2 the read's
    segment at the CR (for non-crossers: of its alignments overlapping it)."""
    near = [a for a in alns if a["chr"] == chrom and a["re"] > s - NEAR and a["rs"] < e + NEAR]
    single = [a for a in near if a["rs"] <= s and a["re"] >= e]
    if single:
        a = max(single, key=lambda a: min(s - a["rs"], a["re"] - e))
        return ("single", min(s - a["rs"], a["re"] - e), *read_interval(a, a, s, e))
    left = [a for a in near if a["rs"] < s]
    right = [a for a in near if a["re"] > e]
    best, pair = -1, None
    for a in left:
        for b in right:
            if a is b or a["strand"] != b["strand"]:
                continue
            # collinear on the read: forward a before b, reverse b before a
            ok = a["qe"] <= b["qs"] + 500 if a["strand"] == "+" else b["qe"] <= a["qs"] + 500
            if ok and b["re"] >= a["re"] and min(s - a["rs"], b["re"] - e) > best:
                best, pair = min(s - a["rs"], b["re"] - e), (a, b)
    if best >= 0:
        return ("split", best, *read_interval(pair[0], pair[1], s, e))
    here = [a for a in near if a["rs"] < e and a["re"] > s]
    qs = [q for a in here for q in read_interval(a, a, s, e)]
    q1, q2 = (min(qs), max(qs)) if qs else (-1, -1)
    if left and right:
        return "both_sides", -1, q1, q2
    return ("left" if left else "right" if right else "inside"), -1, q1, q2


def study(cid: int, crs: list, cache: RC.ReadSequenceCache, rng: random.Random,
          sample: int) -> tuple[dict, list[dict], dict[str, str]]:
    reads, seqs, n_alns, n_supp, alen = [], {}, 0, 0, []
    for cr in crs:
        cr_alns, cr_seqs = cache.fetch_for_cr(cr)
        s, e = cr.referenceStart, cr.referenceEnd
        n_alns += len(cr_alns)
        n_supp += sum(a.is_supplementary for a in cr_alns)
        alen += [a.query_alignment_length for a in cr_alns]
        by_read: dict[str, list] = {}
        for a in cr_alns:
            by_read.setdefault(a.query_name, []).append(a)
        rows = []
        for rn, al in by_read.items():
            prim = [a for a in al if not a.is_supplementary] or al
            kind, anchor, q1, q2 = crossing(sa_alignments(prim[0]), cr.chr, s, e)
            local = max(a.query_alignment_length for a in al)
            rows.append({"container": cid, "crID": cr.crID, "read": rn, "kind": kind,
                         "anchor": anchor, "q1": q1, "q2": q2, "n_alns_here": len(al), "longest_aln_here": local,
                         "read_len": len(cr_seqs[rn].seq) if rn in cr_seqs else 0})
        crossers = [r for r in rows if r["anchor"] >= 0]
        others = [r for r in rows if r["anchor"] < 0]
        picked = crossers + rng.sample(others, min(sample, len(others)))
        for r in rows:
            r["mapped"] = r in picked
        for r in picked:
            if r["read"] in cr_seqs:
                seqs[r["read"]] = str(cr_seqs[r["read"]].seq)
        reads += rows
    rd = pd.DataFrame(reads)
    cross = rd[rd.anchor >= 0].drop_duplicates("read")
    frag = rd[(rd.anchor < 0) & (rd.longest_aln_here < 2000)]
    feats = {
        "container": cid, "chr": crs[0].chr, "start": min(c.referenceStart for c in crs),
        "span": sum(c.referenceEnd - c.referenceStart for c in crs), "n_crs": len(crs),
        "alns": n_alns, "reads": rd.read.nunique(), "supp_share": n_supp / max(1, n_alns),
        "alns_per_read": n_alns / max(1, rd.read.nunique()),
        "median_aln_len": float(np.median(alen)) if alen else 0.0,
        "single": int((cross.kind == "single").sum()), "split": int((cross.kind == "split").sum()),
        "anchor_median": float(cross.anchor.median()) if len(cross) else 0.0,
        "frag_reads": frag.read.nunique(), "frag_share": frag.read.nunique() / max(1, rd.read.nunique()),
    }
    return feats, reads, seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("variant")
    ap.add_argument("out", type=Path)
    ap.add_argument("--min-seconds", type=float, default=120)
    ap.add_argument("--mid-seconds", type=float, default=60)
    ap.add_argument("--controls", type=int, default=100)
    ap.add_argument("--sample", type=int, default=60)
    ap.add_argument("--extra", default="")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(1)
    base = f"{a.root}/results/{a.variant}/20x/HG002"
    c = CV.containers(f"{base}/work")
    c["repr"] = CV.representation(f"{base}/consensus_q100.tsv").reindex(c.index)
    slow = c[c.seconds > a.min_seconds].index.tolist()
    mid = c[(c.seconds > a.mid_seconds) & (c.seconds <= a.min_seconds)].index.tolist()
    fast = c[(c.seconds < 20) & (c.repr >= 0.99)].index.tolist()
    groups = {**dict.fromkeys(slow, "slow"), **dict.fromkeys(mid, "mid"),
              **dict.fromkeys(rng.sample(fast, min(a.controls, len(fast))), "control")}
    jobs = [(Path(f"{base}/work/crs_containers.db"), cid, g) for cid, g in groups.items()]
    for e in filter(None, a.extra.split(",")):
        db, cid = e.rsplit(":", 1)
        jobs.append((Path(db), int(cid), "extra"))
    feats, reads = [], []
    with open(a.out / "reads.fa", "w") as fa:
        for i, (db, cid, g) in enumerate(jobs):
            crs = C.load_crs_containers_from_db(path_db=db, crIDs=[cid])[cid]["crs"]
            cache = RC.ReadSequenceCache(path_alignments=Path(BAM))
            try:
                f, r, seqs = study(cid, crs, cache, rng, a.sample)
            finally:
                cache.close()
            f["group"] = g
            f["db"] = str(db)
            if g != "extra":
                f["seconds"], f["repr"] = c.seconds.get(cid), c.repr.get(cid)
                f["last_level"] = bool(c["last"].get(cid))
            feats.append(f)
            for x in r:
                x["db"] = str(db)
            reads += r
            for rn, sq in seqs.items():
                fa.write(f">{cid}|{rn}\n{sq}\n")
            print(f"[{i + 1}/{len(jobs)}] {g} {cid}: {f['alns']} alns, {f['reads']} reads, "
                  f"{f['single']} single + {f['split']} split crossers", file=sys.stderr, flush=True)
    pd.DataFrame(feats).to_csv(a.out / "containers.tsv", sep="\t", index=False)
    pd.DataFrame(reads).to_csv(a.out / "reads.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
