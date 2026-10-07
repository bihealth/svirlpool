"""Re-form CRs with svirlpool's own code (signalstrength_to_crs) under
alternative settings and score them. usage: cr_reform.py <variant> [...]"""

import os
import sys
import tempfile
from pathlib import Path

import crtools as T
import pandas as pd

from svirlpool.candidateregions import signalstrength_to_crs as S2C

REF = Path("/home/mayv_c/biodata/local/references/GRCh38/GRCh38.fa")
# formed CR sets are cached here (delete a .crs.db to re-form it)
OUT = Path(os.environ.get("CR_STUDY_OUT", "cr_study_out"))
OUT.mkdir(parents=True, exist_ok=True)
EMPTY_TRF = OUT / "empty.bed"
EMPTY_TRF.write_text("")

VARIANTS = {
    "b1200": {"buffer": 1200, "trf": T.TRF},
    "b600": {"buffer": 600, "trf": T.TRF},
    "b300": {"buffer": 300, "trf": T.TRF},
    "b1200_notrf": {"buffer": 1200, "trf": str(EMPTY_TRF)},
    "b300_notrf": {"buffer": 300, "trf": str(EMPTY_TRF)},
}


def form(name: str, signals: Path = Path(T.W + "signalstrength.tsv.gz"), **kw) -> pd.DataFrame:
    v = VARIANTS.get(name, {}) | kw
    db = OUT / f"{name}.crs.db"
    if not db.exists():
        with tempfile.TemporaryDirectory(dir=OUT) as td:
            S2C.create_candidate_regions(
                reference=REF, signalstrengths=signals, output=db,
                tandem_repeats=Path(v["trf"]), threads=8,
                buffer_region_radius=v["buffer"], min_cr_size=500,
                filter_absolute=v.get("fa", 1.0), filter_normalized=v.get("fn", 0.06),
                cutoff_median_readcount_per_region=6.0, tmp_dir_path=Path(td),
            )
    return T.load_crs(db)


if __name__ == "__main__":
    import logging
    logging.disable(logging.INFO)
    truth = T.load_truth()
    rows = [T.evaluate("current crs.db", T.load_crs(), truth)]
    for name in sys.argv[1:]:
        rows.append(T.evaluate(name, form(name), truth))
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
