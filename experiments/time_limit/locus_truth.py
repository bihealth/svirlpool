"""Step 2 of the read-selection study: which reads of a slow locus carry a
correct allele, judged against the HG002 Q100 v1.1 assembly.

1. True locus, per CR and haplotype: the nearest UNIQUE anchors on either
   side. 5 kbp reference chunks at growing distances from the CR (DISTANCES)
   are mapped to the assembly (minimap2 -x asm20 -c -N 5); a chunk anchors on
   haplotype h when its best hit is on <chr>_h, covers >= 80% of it, and no
   other hit outside the two haplotype copies of <chr> scores >= 90% of it.
   The window is the assembly between the nearest left and right anchors
   (same contig and strand). Satellite flanks are skipped this way, so a
   window can be much longer than the CR (`win_len`).
2. The full reads of locus_reads.py (reads.fa) are mapped to the assembly
   (minimap2 -x map-ont -c --eqx, prebuilt asm20 index).
3. A read is judged on its SEGMENT AT THE CR (locus_reads.py q1-q2, from its
   GRCh38 alignment), the part the consensus assembles:
     resolves  one alignment covers the segment contiguously (+- 100 bp) and
               lies in a haplotype window of the CR; `identity` = its
               identity over the segment (= / (= + X + I + D))
     placed    some alignment of the segment lies in a window
     elsewhere the segment is covered by an alignment outside the windows

Per CR: the crossing reads ranked by anchor, the top k (k = --k-factor x the
controls' median read count) against the other crossers and the sampled
non-crossing reads.

usage: locus_truth.py <locus dir> [--threads 16] [--k-factor 2]
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

ASM = "/home/mayv_c/biodata/local/assemblies/hg002v1.1"
REF = "/home/mayv_c/biodata/local/references/GRCh38/GRCh38.fa"
CHUNK = 5000
DISTANCES = (0, 5_000, 10_000, 20_000, 30_000, 50_000, 75_000, 100_000, 150_000,
             200_000, 300_000)
HAPS = ("MATERNAL", "PATERNAL")
CIG = re.compile(r"(\d+)([=XIDM])")


def minimap(preset: str, query: Path, paf: Path, threads: int, extra=()) -> None:
    if paf.exists():
        return
    tmp = paf.with_suffix(".tmp")
    with open(tmp, "w") as f:
        subprocess.run(["nice", "-n", "10", "minimap2", "-x", preset, "-c", "-t", str(threads),
                        *extra, ASM + ".asm20.mmi", str(query)],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    tmp.rename(paf)


def paf_rows(path: Path):
    with open(path) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            tags = {t[:2]: t[5:] for t in x[12:]}
            yield {"q": x[0], "ql": int(x[1]), "qs": int(x[2]), "qe": int(x[3]),
                   "strand": x[4], "t": x[5], "ts": int(x[7]), "te": int(x[8]),
                   "as": int(tags.get("AS", 0)), "cg": tags.get("cg", "")}


def cr_table(d: Path) -> pd.DataFrame:
    rows = []
    c = pd.read_csv(d / "containers.tsv", sep="\t")
    for db, g in c.groupby("db"):
        con = sqlite3.connect(db)
        for cid in g.container:
            data = con.execute("select data from containers where crID = ?", (int(cid),)).fetchone()[0]
            for cr in json.loads(data)["crs"]:
                rows.append({"run": db.split("/")[-5], "container": cid, "crID": cr["crID"],
                             "chr": cr["chr"], "s": cr["referenceStart"], "e": cr["referenceEnd"]})
    return pd.DataFrame(rows)


def anchor_windows(d: Path, crs: pd.DataFrame, threads: int) -> pd.DataFrame:
    ref = pysam.FastaFile(REF)
    fa = d / "anchors.fa"
    with open(fa, "w") as f:
        for r in crs.itertuples():
            L = ref.get_reference_length(r.chr)
            for dist in DISTANCES:
                for side, a in (("L", r.s - dist - CHUNK), ("R", r.e + dist)):
                    if a < 0 or a + CHUNK > L:
                        continue
                    seq = ref.fetch(r.chr, a, a + CHUNK)
                    if seq.upper().count("N") < CHUNK // 10:
                        f.write(f">{r.run}|{r.crID}|{side}|{dist}\n{seq}\n")
    paf = d / "anchors.paf"
    minimap("asm20", fa, paf, threads, ("--secondary=yes", "-N", "5"))
    hits: dict[str, list] = {}
    for h in paf_rows(paf):
        hits.setdefault(h["q"], []).append(h)
    chrom = {(r.run, r.crID): r.chr for r in crs.itertuples()}
    best: dict[tuple, dict] = {}  # (run, crID, side, hap) -> nearest unique anchor hit
    for q, hs in hits.items():
        run, crID, side, dist = q.split("|")
        key0 = (run, int(crID))
        c = chrom[key0]
        for hap in HAPS:
            own = [h for h in hs if h["t"] == f"{c}_{hap}"]
            if not own:
                continue
            top = max(own, key=lambda h: h["as"])
            if top["qe"] - top["qs"] < 0.8 * top["ql"]:
                continue
            rivals = [h for h in hs if h is not top and not h["t"].startswith(f"{c}_")]
            rivals += [h for h in own if h is not top]
            if any(h["as"] >= 0.9 * top["as"] for h in rivals):
                continue
            key = (*key0, side, hap)
            if key not in best or int(dist) < best[key]["dist"]:
                best[key] = {**top, "dist": int(dist)}
    rows = []
    for r in crs.itertuples():
        for hap in HAPS:
            a, b = best.get((r.run, r.crID, "L", hap)), best.get((r.run, r.crID, "R", hap))
            if not a or not b or a["t"] != b["t"] or a["strand"] != b["strand"]:
                continue
            w0, w1 = (a["te"], b["ts"]) if a["strand"] == "+" else (b["te"], a["ts"])
            if not 0 < w1 - w0 < 2_000_000:
                continue
            rows.append({"run": r.run, "crID": r.crID, "hap": hap, "t": a["t"], "w0": w0,
                         "w1": w1, "win_len": w1 - w0, "cr_len": r.e - r.s,
                         "dist_L": a["dist"], "dist_R": b["dist"]})
    return pd.DataFrame(rows)


def segment_identity(h: dict, q1: int, q2: int) -> float:
    """Identity of a PAF alignment (--eqx) over the query interval [q1, q2)."""
    q = h["qs"] if h["strand"] == "+" else h["qe"]
    step = 1 if h["strand"] == "+" else -1
    m = n = 0
    for k, op in CIG.findall(h["cg"]):
        k = int(k)
        if op in "=XMI":
            lo, hi = sorted((q, q + step * k))
            ov = max(0, min(hi, q2) - max(lo, q1))
            if op == "=":
                m += ov
            n += ov
            q += step * k
        elif op == "D" and q1 <= q < q2:
            n += k
    return m / n if n else np.nan


def judge(paf: Path, reads: pd.DataFrame, win: pd.DataFrame, chrom: dict) -> pd.DataFrame:
    """`chrom_ok` (no window needed): one alignment on a haplotype contig of
    the CR's chromosome covers the segment contiguously; `chrom_identity`."""
    segs: dict[str, list] = {}
    for r in reads.itertuples():
        if r.q1 >= 0 and r.q2 > r.q1:
            segs.setdefault(f"{r.container}|{r.read}", []).append((r.run, r.crID, r.q1, r.q2))
    wins: dict[tuple, list] = {}
    for w in win.itertuples():
        wins.setdefault((w.run, w.crID), []).append(w)
    out: dict[tuple, dict] = {}
    for h in paf_rows(paf):
        for run, crID, q1, q2 in segs.get(h["q"], ()):
            if h["qe"] <= q1 or h["qs"] >= q2:
                continue
            key = (run, crID, h["q"].split("|", 1)[1])
            o = out.setdefault(key, {"placed": False, "resolves": False, "elsewhere": False,
                                     "hap": None, "identity": np.nan, "chrom_ok": False,
                                     "chrom_hap": None, "chrom_identity": np.nan})
            covers = h["qs"] <= q1 + 100 and h["qe"] >= q2 - 100
            if covers and h["t"].startswith(f"{chrom[(run, crID)]}_"):
                ident = segment_identity(h, q1, q2)
                if not o["chrom_ok"] or ident > o["chrom_identity"]:
                    o.update(chrom_ok=True, chrom_hap=h["t"].rsplit("_", 1)[1],
                             chrom_identity=ident)
            inside = [w for w in wins.get((run, crID), ())
                      if w.t == h["t"] and h["ts"] < w.w1 and h["te"] > w.w0]
            if not inside:
                o["elsewhere"] = True
                continue
            o["placed"] = True
            if h["qs"] <= q1 + 100 and h["qe"] >= q2 - 100:
                ident = segment_identity(h, q1, q2)
                if not o["resolves"] or ident > o["identity"]:
                    o.update(resolves=True, hap=inside[0].hap, identity=ident)
    return pd.DataFrame([{"run": k[0], "crID": k[1], "read": k[2], **v} for k, v in out.items()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", type=Path)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--k-factor", type=float, default=2.0)
    a = ap.parse_args()
    d = a.dir
    cont = pd.read_csv(d / "containers.tsv", sep="\t")
    reads = pd.read_csv(d / "reads.tsv", sep="\t")
    reads["run"] = reads.db.str.split("/").str[-5]
    crs = cr_table(d)
    win = anchor_windows(d, crs, a.threads)
    win.to_csv(d / "windows.tsv", sep="\t", index=False)
    rpaf = d / "reads.paf"
    minimap("map-ont", d / "reads.fa", rpaf, a.threads, ("--eqx",))
    mapped = reads[reads.mapped]
    j = judge(rpaf, mapped, win, {(x.run, x.crID): x.chr for x in crs.itertuples()})
    r = mapped.merge(j, on=["run", "crID", "read"], how="left")
    for col in ("placed", "resolves", "elsewhere", "chrom_ok"):
        r[col] = r[col].astype("boolean").fillna(False).astype(bool)
    ctrl = cont[cont.group == "control"].reads.median()
    k = int(round(a.k_factor * ctrl))
    r["rank"] = r[r.anchor >= 0].groupby(["run", "crID"]).anchor.rank(ascending=False, method="first")
    r["set"] = np.where(r.anchor < 0, "non-crossing",
                        np.where(r["rank"] <= k, f"top{k} crossing", "other crossing"))
    r.to_csv(d / "reads_judged.tsv", sep="\t", index=False)
    nw = win.groupby(["run", "crID"]).size().rename("n_wins")
    meta = cont.set_index("container")
    crs["group"] = crs.container.map(meta.group)
    crs = crs.join(nw, on=["run", "crID"])
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.3f}"  # noqa: E731
    print(f"k = {k} ({a.k_factor} x the controls' median read count {ctrl:.0f})")
    print("\n== CRs with a window on 0 / 1 / 2 haplotypes; window length / CR length")
    crs["n_wins"] = crs.n_wins.fillna(0).astype(int)
    w = win.assign(group=win.crID.map(crs.set_index("crID").group))
    print(pd.concat([crs.groupby("group").n_wins.value_counts().unstack(fill_value=0),
                     w.groupby("group").apply(lambda g: (g.win_len / g.cr_len).median())
                     .rename("median win/cr")], axis=1).to_string(float_format=fmt))
    r = r.join(nw, on=["run", "crID"])
    r["group"] = r.container.map(meta.group)
    print("\n== ALL CRs, chromosome-level check (no window): segment contiguous on a "
          "haplotype contig of the CR's chromosome")
    print(r.groupby(["group", "set"]).agg(reads=("read", "size"), chrom_ok=("chrom_ok", "mean"),
                                          identity=("chrom_identity", "median"))
          .to_string(float_format=fmt))
    rows = []
    for _, x in r.groupby(["run", "crID"]):
        top = x[x.set.str.startswith("top")]
        rest = x[~x.set.str.startswith("top")]
        rows.append({"group": x.group.iloc[0], "crossers": int((x.anchor >= 0).sum()),
                     "top_ok_share": top.chrom_ok.mean() if len(top) else np.nan,
                     "rest_ok_share": rest.chrom_ok.mean() if len(rest) else np.nan,
                     "haps_all": x[x.chrom_ok].chrom_hap.nunique(),
                     "haps_top": top[top.chrom_ok].chrom_hap.nunique(),
                     "no_crossers": len(top) == 0})
    q = pd.DataFrame(rows)
    print(q.groupby("group").apply(lambda g: pd.Series({
        "CRs": len(g), "no_crossers": g.no_crossers.mean(),
        "top_ok_share": g.top_ok_share.mean(), "rest_ok_share": g.rest_ok_share.mean(),
        "both_haps_all": (g.haps_all == 2).mean(), "both_haps_top": (g.haps_top == 2).mean()}))
        .to_string(float_format=fmt))
    r = r[r.n_wins.notna()]
    g = r.groupby(["group", "set"]).agg(reads=("read", "size"), placed=("placed", "mean"),
                                        resolves=("resolves", "mean"),
                                        elsewhere=("elsewhere", "mean"),
                                        identity=("identity", "median"))
    print("\n== reads of CRs with a window: share placed / resolving / elsewhere")
    print(g.to_string(float_format=fmt))
    rows = []
    for (run, crID), x in r.groupby(["run", "crID"]):
        top = x[x.set.str.startswith("top")]
        rest = x[~x.set.str.startswith("top")]
        rows.append({"run": run, "crID": crID, "container": x.container.iloc[0],
                     "group": x.group.iloc[0], "n_wins": int(x.n_wins.iloc[0]),
                     "crossers": int((x.anchor >= 0).sum()), "top": len(top),
                     "top_res": int(top.resolves.sum()), "rest_res": int(rest.resolves.sum()),
                     "haps_all": x[x.resolves].hap.nunique(),
                     "haps_top": top[top.resolves].hap.nunique(),
                     "top_res_share": top.resolves.mean() if len(top) else np.nan,
                     "rest_res_share": rest.resolves.mean() if len(rest) else np.nan})
    p = pd.DataFrame(rows)
    p.to_csv(d / "crs_judged.tsv", sep="\t", index=False)
    print("\n== per CR with windows on both haplotypes: resolving reads, haplotypes covered")
    p2 = p[p.n_wins == 2]
    print(p2.groupby("group")[["crossers", "top", "top_res", "rest_res", "haps_all", "haps_top",
                               "top_res_share", "rest_res_share"]].mean()
          .to_string(float_format=fmt))
    print(p2.groupby("group").apply(lambda g: pd.Series({
        "CRs": len(g), "both_haps_all": (g.haps_all == 2).mean(),
        "both_haps_top": (g.haps_top == 2).mean(),
        "top_lost_a_hap": (g.haps_top < g.haps_all).mean()})).to_string(float_format=fmt))


if __name__ == "__main__":
    main()
