"""Similarity of every core consensus of a svirlpool run to its counterpart in
the HG002 Q100 v1.1 diploid assembly (the T2TQ100 truth), after the H1
measurement of svp_merging (workflow/scripts/h1/measure.py):

1. the PADDED consensus (consensus.fasta: the core plus raw-read padding) is
   mapped to the assembly (minimap2 -x asm20) only to find the locus;
2. on each haplotype (MATERNAL, PATERNAL; contigs of the consensus' reference
   chromosome preferred), the core's two ends are projected through the
   best chain of alignments, and the assembly sequence between them, plus a
   margin, is the target;
3. the CORE ALONE is aligned to it with edlib in infix mode (the whole core
   aligns, the target's ends are free). No padding base is scored.

Per consensus: edit distance to either haplotype (ed_MATERNAL, ed_PATERNAL),
the closer one (hap, ed) and identity = 1 - ed / len(core).

usage: consensus_q100.py <workdir> <out.tsv> [--threads 24]
       [--fasta hg002v1.1.fasta] [--mmi hg002v1.1.asm20.mmi]
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import shutil
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path

import edlib
import pandas as pd
import pysam

ASM = Path("/home/mayv_c/biodata/local/assemblies")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTN", "TGCAN")
HAPS = ("MATERNAL", "PATERNAL")
FASTA = None


def revcomp(s):
    return s.translate(COMP)[::-1]


def load(work: Path) -> dict:
    out = {}
    for f in sorted(glob.glob(str(work / "consensus" / "*" / "*.jsonl"))):
        for line in open(f):
            o = json.loads(line)
            for cid, c in o["consensus_dicts"].items():
                pad = c["consensus_padding"]
                core = c["consensus_sequence"].upper()
                L = pad["padding_size_left"]
                regs = c.get("original_regions") or []
                chrom = (
                    regs[0][0] if regs and isinstance(regs[0], (list, tuple)) else None
                )
                meta = c.get("clustering_meta_data") or {}
                out[cid] = {
                    "core": core,
                    "padded": pad["sequence"],
                    "L": L,
                    "chrom": chrom,
                    "phasing_status": meta.get("phasing_status"),
                    "phasing_sites": meta.get("phasing_sites"),
                    "n_alleles": meta.get("n_alleles"),
                }
    return out


def read_paf(path):
    by_q = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tags = dict(t.split(":", 2)[::2] for t in f[12:])
            by_q.setdefault(f[0], []).append(
                {
                    "qs": int(f[2]),
                    "qe": int(f[3]),
                    "strand": f[4],
                    "t": f[5],
                    "ts": int(f[7]),
                    "te": int(f[8]),
                    "as": int(tags.get("AS", 0)),
                    "cg": tags.get("cg", ""),
                }
            )
    return by_q


def project(aln, q):
    """Target coordinate of query position q through one PAF alignment, or None."""
    if not aln["qs"] <= q <= aln["qe"]:
        return None
    t = aln["ts"]
    qpos = aln["qs"] if aln["strand"] == "+" else aln["qe"]
    step = 1 if aln["strand"] == "+" else -1
    for n, op in CIGAR_RE.findall(aln["cg"]):
        n = int(n)
        if op in "M=X":
            qnext = qpos + step * n
            if min(qpos, qnext) <= q <= max(qpos, qnext):
                return t + abs(q - qpos)
            qpos, t = qnext, t + n
        elif op == "I":
            qnext = qpos + step * n
            if min(qpos, qnext) <= q <= max(qpos, qnext):
                return t
            qpos = qnext
        elif op in "DN":
            t += n
    return t


def locus(alns, hap, chrom, q0, q1):
    cand = [a for a in alns if a["t"].endswith(hap)]
    same = [a for a in cand if chrom and a["t"] == f"{chrom}_{hap}"]
    cand = same or cand
    if not cand:
        return None
    groups = {}
    for a in cand:
        groups.setdefault((a["t"], a["strand"]), []).append(a)
    best, best_score = None, -1
    for key, g in groups.items():
        top = max(g, key=lambda a: a["as"])
        chain = [a for a in g if abs(a["ts"] - top["ts"]) < 300_000]
        score = sum(
            a["as"] for a in chain if a["qe"] > q0 - 5000 and a["qs"] < q1 + 5000
        )
        if score > best_score:
            best, best_score = (key, chain), score
    (tname, strand), chain = best
    ends, extrap = [], False
    for q in (q0, q1):
        hits = [p for a in chain if (p := project(a, q)) is not None]
        if hits:
            ends.append(hits[0])
            continue
        extrap = True
        a = min(chain, key=lambda a: min(abs(q - a["qs"]), abs(q - a["qe"])))
        ends.append(
            a["ts"] + (q - a["qs"]) if strand == "+" else a["te"] - (q - a["qs"])
        )
    return tname, strand, min(ends), max(ends), extrap


def init(fasta):
    global FASTA
    FASTA = pysam.FastaFile(fasta)


def one(job):
    cid, core, L, chrom, alns = job
    res = {"consensus": cid, "len": len(core), "status": "ok"}
    if not alns:
        res["status"] = "unmapped"
        return res
    q0, q1 = L, L + len(core)
    for hap in HAPS:
        loc = locus(alns, hap, chrom, q0, q1)
        if loc is None:
            continue
        tname, strand, t0, t1, extrap = loc
        margin = 500 + (len(core) if extrap else len(core) // 5)
        tlen = FASTA.get_reference_length(tname)
        t0, t1 = min(max(0, t0 - margin), tlen), min(max(0, t1 + margin), tlen)
        if t1 - t0 < len(core) // 2:
            continue
        target = FASTA.fetch(tname, t0, t1).upper()
        if strand == "-":
            target = revcomp(target)
        r = edlib.align(core, target, mode="HW", task="distance")
        res[f"ed_{hap}"] = r["editDistance"]
        res[f"t_{hap}"] = f"{tname}:{t0}-{t1}{strand}"
        res[f"extrap_{hap}"] = extrap
    eds = {h: res[f"ed_{h}"] for h in HAPS if f"ed_{h}" in res}
    if not eds:
        res["status"] = "no_locus"
        return res
    hap = min(eds, key=eds.get)
    res["hap"], res["ed"] = hap, eds[hap]
    res["identity"] = 1 - eds[hap] / max(1, len(core))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--fasta", default=str(ASM / "hg002v1.1.fasta"))
    ap.add_argument("--mmi", default=str(ASM / "hg002v1.1.asm20.mmi"))
    a = ap.parse_args()
    cons = load(a.workdir)
    tmp = Path(tempfile.mkdtemp(prefix="cq100_", dir=a.out.parent))
    try:
        fa = tmp / "padded.fa"
        with open(fa, "w") as f:
            for cid, c in cons.items():
                f.write(f">{cid}\n{c['padded']}\n")
        paf = tmp / "padded_vs_q100.paf"
        with open(paf, "w") as f:
            subprocess.run(
                [
                    "minimap2",
                    "-c",
                    "-t",
                    str(a.threads),
                    "--secondary=yes",
                    "-N",
                    "5",
                    a.mmi,
                    str(fa),
                ],
                stdout=f,
                stderr=subprocess.DEVNULL,
                check=True,
            )
        alns = read_paf(paf)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    jobs = [
        (cid, c["core"], c["L"], c["chrom"], alns.get(cid, []))
        for cid, c in cons.items()
    ]
    jobs.sort(key=lambda j: -len(j[1]))  # long cores first: better load balance
    with Pool(a.threads, initializer=init, initargs=(a.fasta,)) as pool:
        rows = pool.map(one, jobs, chunksize=4)
    d = pd.DataFrame(rows)
    meta = pd.DataFrame(
        [
            {
                "consensus": k,
                "chrom": c["chrom"],
                "phasing_status": c["phasing_status"],
                "phasing_sites": c["phasing_sites"],
                "n_alleles": c["n_alleles"],
            }
            for k, c in cons.items()
        ]
    )
    d = meta.merge(d, on="consensus")
    d["crID"] = d.consensus.str.split(".").str[0].astype(int)
    d.sort_values(["crID", "consensus"]).to_csv(a.out, sep="\t", index=False)
    ok = d[d.status == "ok"]
    print(
        f"{len(d)} consensuses, {d.status.value_counts().to_dict()}; "
        f"identity mean {ok.identity.mean():.4f}, aggregate {1 - ok.ed.sum() / ok.len.sum():.4f}"
    )


if __name__ == "__main__":
    main()
