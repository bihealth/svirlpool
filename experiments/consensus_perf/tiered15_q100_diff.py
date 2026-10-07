"""Where does ps_ref beat ps_ava against Q100? Per container: ava's phasing
status and consensus count vs ref's, and the repr / prec differences."""

import sys

import pandas as pd

sys.path.insert(
    0,
    "/home/mayv_c/development/svirlpool/.claude/worktrees/ref-vs-ava/experiments/consensus_perf",
)
from tiered15_report import representation  # noqa: E402

R = "/home/mayv_c/development/svp_tiered15/results"
a = pd.read_csv(f"{R}/ps_ava/20x/HG002/consensus_q100.tsv", sep="\t")
b = pd.read_csv(
    f"{R}/{sys.argv[1] if len(sys.argv) > 1 else 'ps_ref'}/20x/HG002/consensus_q100.tsv",
    sep="\t",
)
ra, rb = representation(a), representation(b)


def st(d):
    return d.groupby("crID").agg(
        status=("phasing_status", "first"),
        n=("consensus", "size"),
        n_alleles=("n_alleles", "first"),
        minlen=("len", "min"),
    )


j = (
    ra[["repr", "prec"]]
    .join(st(a))
    .join(rb[["repr", "prec"]].join(st(b)), rsuffix="_b", how="inner")
)
j["d_repr"], j["d_prec"] = j.repr_b - j.repr, j.prec_b - j.prec
j["cls"] = (
    j.status.astype(str)
    + " "
    + j.n.astype(str)
    + " -> "
    + j.status_b.astype(str)
    + " "
    + j.n_b.astype(str)
)
diff = j[(j.d_repr.abs() > 0.005) | (j.d_prec.abs() > 0.005)]
pd.set_option("display.width", 200)
print(f"{len(diff)} of {len(j)} containers differ by > 0.005")
g = diff.groupby("cls").agg(
    n=("d_repr", "size"),
    d_repr=("d_repr", "sum"),
    d_prec=("d_prec", "sum"),
    repr_better=("d_repr", lambda x: (x > 0.005).sum()),
    repr_worse=("d_repr", lambda x: (x < -0.005).sum()),
)
print(
    g.sort_values("n", ascending=False)
    .head(15)
    .to_string(float_format=lambda x: f"{x:.3f}")
)
print("\ncontribution to the mean repr difference:", round(j.d_repr.sum() / len(j), 5))
print("\nlargest ref gains:")
print(
    diff.sort_values("d_repr")
    .tail(8)[["cls", "repr", "repr_b", "prec", "prec_b", "minlen"]]
    .to_string()
)

# by consensus batch (genomic block)
bt = pd.read_csv(f"{R}/ps_ava/20x/HG002/work/consensus_batches.tsv", sep="\t")
loc = {
    int(c): (r.chr, r.start, r.end)
    for r in bt.itertuples()
    for c in str(r.crIDs).split(",")
}
j["chr"] = [loc[c][0] if c in loc else "?" for c in j.index]
j["mb"] = [loc[c][1] // 1_000_000 if c in loc else -1 for c in j.index]
j["block"] = [
    f"{loc[c][0]}:{loc[c][1] // 1_000_000}-{loc[c][2] // 1_000_000}Mb"
    if c in loc
    else "?"
    for c in j.index
]
pb = j.groupby("block").agg(
    n=("d_repr", "size"), d_repr=("d_repr", "sum"), d_prec=("d_prec", "sum")
)
pb["share"] = pb.d_repr / j.d_repr.sum()
print("\nbatches by contribution to the repr difference:")
print(
    pb.sort_values("d_repr", ascending=False)
    .head(8)
    .to_string(float_format=lambda x: f"{x:.3f}")
)
j["end_mb"] = [loc[c][2] // 1_000_000 if c in loc else -1 for c in j.index]
cen = (j.chr == "chr19") & (j.end_mb >= 24) & (j.mb < 30)
x = j[~cen]
print(
    f"\nwithout chr19 24-30 Mb ({cen.sum()} containers): repr diff {x.d_repr.mean():.5f}, "
    f"prec diff {x.d_prec.mean():.5f}, better/worse repr "
    f"{(x.d_repr > 0.005).sum()}/{(x.d_repr < -0.005).sum()}, prec "
    f"{(x.d_prec > 0.005).sum()}/{(x.d_prec < -0.005).sum()}"
)
