"""End-to-end comparison of --phasing-sites on the trio over the
representative 15% region set (svp_tiered15, a copy of svp_improvements with
config/regions15.bed and truvari --pctsize 0.9).

Per variant:
  sv        truvari refined P / R / F1 (V5, T2TQ100; all, non_trf), raw F1
  mendel    HG002 trio Mendelian consistency (all SVs >= 50 bp)
  time      svirlpool run per sample: wall clock, CPU (user + sys of the whole
            run, children included), peak RSS, summed consensus batch seconds,
            containers escalated after a tool timeout and consensuses that
            still timed out at the last level (wall-clock timeouts depend on
            the machine's load: a confounder between variants)
  q100      HG002 core consensuses vs the Q100 assembly (consensus_q100.py
            output <results>/<variant>/20x/HG002/consensus_q100.tsv):
            located share, mean and aggregate identity, share >= 0.99 / 0.999
Paired by container against the reference variant (--ref, default ps_ava):
  per container and Q100 haplotype, the best identity of any of its
  consensuses (1 - ed_h / len); a container's representation identity is the
  mean over the two haplotypes. Containers better / worse by > 0.005 and the
  mean difference, overall and, for the tiered variant, split by the sites
  its phasing used (reference / ava).

usage: tiered15_report.py <svp_tiered15 dir> [--variants ps_ava,ps_tiered,ps_ref] [--ref ps_ava]
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import pandas as pd

HAPS = ("MATERNAL", "PATERNAL")
SAMPLES = ("HG002", "HG003", "HG004")


def sv_table(root, variants):
    rows = []
    for v in variants:
        for ts in ("V5", "T2TQ100"):
            for rs in ("all", "non_trf"):
                d = root / "results" / v / "20x" / "truvari" / ts / rs
                ref = d / "refine.variant_summary.json"
                raw = d / "summary.json"
                if not ref.exists():
                    continue
                r = json.load(open(ref))
                w = json.load(open(raw)) if raw.exists() else {}
                rows.append(
                    {
                        "variant": v,
                        "truthset": ts,
                        "regions": rs,
                        "P": r["precision"],
                        "R": r["recall"],
                        "F1": r["f1"],
                        "TP": r["TP-base"],
                        "FP": r["FP"],
                        "FN": r["FN"],
                        "F1_raw": w.get("f1"),
                        "GT_conc_raw": w.get("gt_concordance"),
                    }
                )
    return pd.DataFrame(rows)


def mendel_table(root, variants):
    rows = []
    for v in variants:
        f = root / "results" / v / "20x" / "mendel" / "mendel.tsv"
        if not f.exists():
            continue
        m = pd.read_csv(f, sep="\t")
        m = m[(m.svtype == "all") & (m.size_bin == "all")]
        n = m.set_index("status")["count"]
        c, i = int(n["consistent"]), int(n["inconsistent"])
        rows.append(
            {"variant": v, "consistent": c, "inconsistent": i, "MC": c / (c + i)}
        )
    return pd.DataFrame(rows)


def _time_txt(path):
    t = open(path).read()
    if "User time" not in t:  # the run is not finished
        return None
    user = float(re.search(r"User time \(seconds\): ([\d.]+)", t).group(1))
    sys_ = float(re.search(r"System time \(seconds\): ([\d.]+)", t).group(1))
    el = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t
    ).group(1)
    secs = 0.0
    for part in el.split(":"):
        secs = secs * 60 + float(part)
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))
    return secs, user + sys_, rss / 1e6


def time_table(root, variants):
    rows = []
    for v in variants:
        for s in SAMPLES:
            f = root / "benchmarks" / "svirlpool" / f"run.{v}.20x.{s}.time.txt"
            if not f.exists():
                continue
            parsed = _time_txt(f)
            if parsed is None:
                continue
            wall, cpu, rss = parsed
            cons = 0.0
            for b in glob.glob(
                str(
                    root
                    / "results"
                    / v
                    / "20x"
                    / s
                    / "work"
                    / "benchmarks"
                    / "consensus"
                    / "*.txt"
                )
            ):
                cons += float(pd.read_csv(b, sep="\t")["s"].iloc[0])
            esc, last = 0, set()
            for lg in glob.glob(
                str(
                    root
                    / "results"
                    / v
                    / "20x"
                    / s
                    / "work"
                    / "consensus"
                    / "*"
                    / "*.log"
                )
            ):
                if lg.endswith(".diag.log"):
                    continue
                for line in open(lg):
                    if m := re.search(
                        r"(\d+) container\(s\) escalated after a tool timeout", line
                    ):
                        esc += int(m.group(1))
                    elif m := re.search(r"timed out for (\S+) after 120 seconds", line):
                        last.add(m.group(1))
            rows.append(
                {
                    "variant": v,
                    "sample": s,
                    "wall_min": wall / 60,
                    "cpu_h": cpu / 3600,
                    "rss_gb": rss,
                    "consensus_batch_h": cons / 3600,
                    "escalated": esc,
                    "timeout_last_level": len(last),
                }
            )
    return pd.DataFrame(rows)


def load_q100(root, v):
    f = root / "results" / v / "20x" / "HG002" / "consensus_q100.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


def q100_summary(d):
    ok = d[d.status == "ok"]
    return {
        "consensuses": len(d),
        "located": len(ok) / max(1, len(d)),
        "identity_mean": ok.identity.mean(),
        "identity_aggr": 1 - ok.ed.sum() / ok.len.sum(),
        "ge_0.99": (ok.identity >= 0.99).mean(),
        "ge_0.999": (ok.identity >= 0.999).mean(),
        "containers": d.crID.nunique(),
        "cons_per_container": len(d) / max(1, d.crID.nunique()),
    }


def representation(d):
    """Per container: best identity to each Q100 haplotype and their mean
    (repr), and the mean identity of its consensuses to their closer
    haplotype (prec)."""
    ok = d[d.status == "ok"].copy()
    for h in HAPS:
        ok[f"id_{h}"] = 1 - ok[f"ed_{h}"] / ok.len.clip(lower=1)
    g = ok.groupby("crID")
    r = pd.DataFrame({f"best_{h}": g[f"id_{h}"].max() for h in HAPS})
    r["repr"] = r[[f"best_{h}" for h in HAPS]].mean(axis=1, skipna=True)
    r["prec"] = g.identity.mean()
    r["n_cons"] = g.size()
    r["sites"] = (
        g.phasing_sites.agg(lambda s: s.dropna().iloc[0] if s.notna().any() else None)
        if "phasing_sites" in ok
        else None
    )
    return r


def timed_out_containers(root, v, sample="HG002"):
    """crIDs with a consensus that timed out at the last escalation level."""
    out = set()
    for lg in glob.glob(
        str(
            root / "results" / v / "20x" / sample / "work" / "consensus" / "*" / "*.log"
        )
    ):
        if lg.endswith(".diag.log"):
            continue
        for line in open(lg):
            if m := re.search(r"timed out for (\d+)\.\d+ after 120 seconds", line):
                out.add(int(m.group(1)))
    return out


def _cmp(j, thr):
    out = {"containers": len(j)}
    for m in ("repr", "prec"):
        d = j[m] - j[f"{m}_ref"]
        out[m] = {
            "diff": round(float(d.mean()), 5),
            "better": int((d > thr).sum()),
            "worse": int((d < -thr).sum()),
        }
    out["cons_per_container"] = (
        round(float(j.n_cons_ref.mean()), 3),
        round(float(j.n_cons.mean()), 3),
    )
    out["n_cons_changed"] = int((j.n_cons != j.n_cons_ref).sum())
    return out


def paired(ref, other, thr=0.005):
    """Containers of both variants: repr (each Q100 haplotype's best match,
    recall-like) and prec (mean identity of the container's consensuses,
    precision-like), better / worse by > thr."""
    j = ref[["repr", "prec", "n_cons"]].join(
        other[["repr", "prec", "n_cons", "sites"]],
        lsuffix="_ref",
        rsuffix="",
        how="inner",
    )
    j = j.dropna(subset=["repr_ref", "repr"])
    groups = {}
    if j.sites.notna().any():
        for site, g in j.groupby("sites"):
            groups[site] = _cmp(g, thr)
    return _cmp(j, thr), groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--variants", default="ps_ava,ps_tiered,ps_ref")
    ap.add_argument("--ref", default="ps_ava")
    a = ap.parse_args()
    variants = a.variants.split(",")
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    fmt = lambda x: f"{x:.4f}"  # noqa: E731
    print("== SV benchmark (truvari --pctsize 0.9, refined)")
    sv = sv_table(a.root, variants)
    if len(sv):
        print(sv.to_string(index=False, float_format=fmt))
    print("\n== Mendelian consistency (trio)")
    print(mendel_table(a.root, variants).to_string(index=False, float_format=fmt))
    print("\n== svirlpool run")
    t = time_table(a.root, variants)
    if len(t):
        print(t.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
        print(
            t.groupby("variant")[
                [
                    "wall_min",
                    "cpu_h",
                    "consensus_batch_h",
                    "escalated",
                    "timeout_last_level",
                ]
            ]
            .sum()
            .to_string(float_format=lambda x: f"{x:.2f}")
        )
    print("\n== HG002 core consensuses vs Q100")
    q = {v: load_q100(a.root, v) for v in variants}
    q = {v: d for v, d in q.items() if d is not None}
    if q:
        print(
            pd.DataFrame({v: q100_summary(d) for v, d in q.items()}).T.to_string(
                float_format=fmt
            )
        )
        for v, d in q.items():
            if "phasing_sites" in d and d.phasing_sites.notna().any():
                c = d.groupby("crID").phasing_sites.first().value_counts()
                print(f"  {v}: containers by phasing sites {c.to_dict()}")
    if a.ref in q:
        rep = {v: representation(d) for v, d in q.items()}
        print(
            f"\n== container representation of the two Q100 haplotypes, paired vs {a.ref}"
        )
        print(
            f"  {a.ref}: mean {rep[a.ref].repr.mean():.5f} over {len(rep[a.ref])} containers"
        )
        for v in q:
            if v == a.ref:
                continue
            o, g = paired(rep[a.ref], rep[v])
            print(f"  {v}: {o}")
            for s, x in g.items():
                print(f"    phased with {s} sites: {x}")
            # wall-clock timeouts depend on the load: without the containers
            # that timed out at the last level in either variant
            drop = timed_out_containers(a.root, a.ref) | timed_out_containers(a.root, v)
            o, g = paired(
                rep[a.ref][~rep[a.ref].index.isin(drop)],
                rep[v][~rep[v].index.isin(drop)],
            )
            print(
                f"  {v}, without {len(drop)} containers with a last-level timeout: {o}"
            )
            for s, x in g.items():
                print(f"    phased with {s} sites: {x}")


if __name__ == "__main__":
    main()
