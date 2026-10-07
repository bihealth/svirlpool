"""Step 2: validate the CR formation rules (--min-signal-support 2, buffer 600,
--satellite-depth-factor) off the tuning set: HG003 / HG004 of the 15% trio run
(no truth: their own trio calls, overall and at Mendel-consistent sites, still
covered by a CR) and HG002 on the disjoint 30% regions (T2TQ100 SVs inside the
benchmark bed). Reports what a variant loses relative to the current CRs."""

import logging
import subprocess
import tempfile
import time
from pathlib import Path

import cr_reform as R
import crtools as T
import numpy as np
import pandas as pd

from svirlpool.candidateregions import signalstrength_to_crs as S2C

logging.disable(logging.INFO)

OUT = R.OUT / "validate"
OUT.mkdir(exist_ok=True)
RUN15 = "/home/mayv_c/development/svp_tiered15/results/ps_tiered/20x/"
RUN30 = "/home/mayv_c/development/phasing_ablation_30pct/"
Q100 = T.S + "resources/truthsets/GRCh38_HG2-T2TQ100-V1.1_stvar"
VARIANTS = {
    "b600": {"buffer": 600},
    "support1 b600": {"buffer": 600, "support": 1},
    "support1 b600 sat1.5": {"buffer": 600, "support": 1, "sat": 1.5},
    "support1 b600 sat1.25": {"buffer": 600, "support": 1, "sat": 1.25},
    "support2 b600": {"buffer": 600, "support": 2},
    "support2 b600 sat1.5": {"buffer": 600, "support": 2, "sat": 1.5},
}


def form(tag: str, work: str, buffer: int, support: int = 0, sat: float = 0.0) -> pd.DataFrame:
    db = OUT / f"{tag}_b{buffer}_s{support}_sat{sat}.crs.db"
    if not db.exists():
        t = time.time()
        with tempfile.TemporaryDirectory(dir=OUT) as td:
            S2C.create_candidate_regions(
                reference=R.REF, signalstrengths=Path(work + "signalstrength.tsv.gz"), output=db,
                tandem_repeats=Path(T.TRF), threads=8, buffer_region_radius=buffer, min_cr_size=500,
                filter_absolute=1.0, filter_normalized=0.06, cutoff_median_readcount_per_region=6.0,
                tmp_dir_path=Path(td), min_signal_support=support, satellite_depth_factor=sat)
        print(f"  formed {db.name} in {time.time() - t:.0f} s", flush=True)
    return T.load_crs(db)


def zcat_rows(path: str, cols: str) -> list[list[str]]:
    out = subprocess.run(["bash", "-c", f"zcat {path} | grep -v '^#' | cut -f{cols}"],
                         capture_output=True, text=True, check=True).stdout
    return [line.split("\t") for line in out.splitlines()]


def alleles(gt: str) -> list[str] | None:
    a = gt.split(":")[0].replace("|", "/").split("/")
    return None if "." in a else a


def trio_calls() -> pd.DataFrame:
    """Per sample: its non-ref calls in family.vcf.gz, with whether the site is
    Mendel-consistent (HG002 child of HG003 x HG004, all three genotyped)."""
    rows = []
    for c, p, ref, alt, _info, *gts in zcat_rows(RUN15 + "family.vcf.gz", "1,2,4,5,8,10,11,12"):
        p = int(p)
        g = [alleles(x) for x in gts]
        size = abs(len(ref) - len(alt))
        end = p + max(1, len(ref) - 1) if len(ref) > len(alt) else p + 1
        mendel = None
        if all(x is not None for x in g) and any(a != "0" for x in g for a in x):
            ch, fa, mo = g
            if len(ch) == 1:  # haploid (sex chromosomes)
                mendel = ch[0] in fa or ch[0] in mo
            else:
                mendel = (ch[0] in fa and ch[1] in mo) or (ch[1] in fa and ch[0] in mo)
        for sample, x in zip(("HG002", "HG003", "HG004"), g, strict=True):
            if x is not None and any(a != "0" for a in x):
                rows.append((sample, c, p, end, size, mendel is True))
    return pd.DataFrame(rows, columns=["sample", "chr", "start", "end", "size", "mendel_ok"])


def truth30() -> pd.DataFrame:
    """T2TQ100 SVs >= 50 bp with both ends in one benchmark interval and in the 30% regions."""
    def bed(path):
        b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        return {c: g.sort_values("s")[["s", "e"]].values for c, g in b.groupby("c")}
    bench, reg = bed(Q100 + ".benchmark.bed"), bed(RUN30 + "regions30.bed")

    def inside(by, c, s, e):
        a = by.get(c)
        if a is None:
            return False
        k = np.searchsorted(a[:, 0], s, side="right") - 1
        return k >= 0 and a[k, 1] >= e

    rows = []
    for c, p, ref, alt in zcat_rows(Q100 + ".vcf.gz", "1,2,4,5"):
        size = abs(len(ref) - len(alt))
        if size < 50:
            continue
        p = int(p)
        end = p + max(1, len(ref) - 1) if len(ref) > len(alt) else p + 1
        if inside(bench, c, p, end) and inside(reg, c, p, end):
            rows.append((c, p, end, size))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "size"])


def cost_row(name: str, crs: pd.DataFrame, heavy: float = 100_000) -> dict:
    c = T.cost(crs)
    span = crs.end - crs.start
    return {"variant": name, "n_crs": len(crs), "n_span_gt5k": int((span > 5000).sum()),
                "cost_Mbp": c.sum() / 1e6, "heavy_cost_Mbp": c[c >= heavy].sum() / 1e6}


def lost(points: pd.DataFrame, base: np.ndarray, crs: pd.DataFrame) -> tuple[int, int]:
    """(points covered by the current CRs, of those no longer covered)"""
    now = T.covered(points, crs)
    return int(base.sum()), int((base & ~now).sum())


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    calls = trio_calls()
    print("trio calls per sample:", calls.groupby("sample").size().to_dict(),
          " at Mendel-consistent sites:", calls[calls.mendel_ok].groupby("sample").size().to_dict())
    rows = []
    for sample in ("HG002", "HG003", "HG004"):
        work = RUN15 + f"{sample}/work/"
        cur = T.load_crs(work + "crs.db")
        own = calls[calls["sample"] == sample].reset_index(drop=True)
        base_all = T.covered(own, cur)
        mok = own.mendel_ok.values
        big = own["size"].values >= 1000
        base_row = cost_row("current", cur)
        rows.append(dict(sample=f"{sample} 15%", **base_row,
                         calls_cov=int(base_all.sum()), calls_lost=0, mendel_ok_lost=0, calls_1kb_lost=0))
        for name, v in VARIANTS.items():
            crs = form(f"{sample}15", work, **v)
            now = T.covered(own, crs)
            gone = base_all & ~now
            rows.append(dict(sample=f"{sample} 15%", **cost_row(name, crs), calls_cov=int(now.sum()),
                             calls_lost=int(gone.sum()), mendel_ok_lost=int((gone & mok).sum()),
                             calls_1kb_lost=int((gone & big).sum())))
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.1f}"))

    truth = truth30()
    work = RUN30 + "HG002/work/"
    cur = T.load_crs(work + "crs.db")
    base = T.covered(truth, cur)
    big = truth["size"].values >= 1000
    rows = [dict(**cost_row("current", cur), truth=len(truth), covered=int(base.sum()), lost=0, lost_1kb=0)]
    check = form("HG002_30", work, buffer=1200)
    assert len(check) == len(cur), (len(check), len(cur))
    for name, v in VARIANTS.items():
        crs = form("HG002_30", work, **v)
        now = T.covered(truth, crs)
        gone = base & ~now
        rows.append(dict(**cost_row(name, crs), truth=len(truth), covered=int(now.sum()),
                         lost=int(gone.sum()), lost_1kb=int((gone & big).sum())))
    print("\nHG002, disjoint 30% regions, T2TQ100 SVs >= 50 bp in the benchmark bed:")
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.1f}"))
