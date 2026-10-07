"""Offline cost of support1 / satellite rule while keeping the 1200 bp merge buffer."""

import cr_validate as V
import crtools as T
import pandas as pd

pd.set_option("display.width", 250)
rows = []
for s in ("HG002", "HG003", "HG004"):
    w = V.RUN15 + f"{s}/work/"
    rows.append(dict(sample=s, **V.cost_row("current b1200", T.load_crs(w + "crs.db"))))
    for name, kw in (("b1200 support1", {"buffer": 1200, "support": 1}),
                     ("b1200 support1 sat1.5", {"buffer": 1200, "support": 1, "sat": 1.5}),
                     ("b600 support1 sat1.5", {"buffer": 600, "support": 1, "sat": 1.5})):
        rows.append(dict(sample=s, **V.cost_row(name, V.form(f"{s}15", w, **kw))))
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.1f}"))
