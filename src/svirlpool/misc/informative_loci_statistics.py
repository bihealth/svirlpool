# python3 ../../../scripts/reads_to_informative_loci.py -i test/g.readAlignmentsSequence.tsv -o test/g.loci.unmerged.bed -r /fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa -m 30 --nomerge
# bedtools intersect -wa -c -a truth.bed -b g.loci.unmerged.bp.bed > g.truth.bp.count.bed
# bedtools intersect -wa -c -a truth.bed -b g.loci.unmerged.in.bed > g.truth.in.count.bed
# paste <(cat truth.bed) <(awk -v OFS='\t' ' { print $5} ' g.truth.bp.count.bed) <(awk -v OFS='\t' ' { print $5} ' g.truth.in.count.bed) > g.truth.counts.bed


# %%
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

# %%

p = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/test/g.truth.counts.bed"
)

df = pd.read_csv(p, sep="\t", header=None)
df.columns = ["chr", "start", "end", "type:size", "breakends", "indels"]
df["SVtype"] = [row.loc["type:size"].split(":")[0] for i, row in df.iterrows()]
df["SVsize"] = [int(row.loc["type:size"].split(":")[1]) for i, row in df.iterrows()]
df["signals"] = df["indels"] + df["breakends"]
df["ii-ratio"] = df["breakends"].div(df["signals"])
df["logSVsize"] = df["SVsize"].apply(np.log)

df.fillna(0, inplace=True)
# %%

df.sort_values(by="SVsize", inplace=True)
fig = px.scatter(
    df,
    x="logSVsize",
    y="signals",
    color="ii-ratio",
    marginal_x="histogram",
    marginal_y="histogram",
    facet_row="SVtype",
    hover_data=["SVsize", "breakends", "indels", "chr", "start"],
)
fig.write_html(p.with_suffix(".html"))

# %%
fig = px.histogram(df, x="SVsize", range_x=(0, 1000), nbins=4000, color="SVtype")
fig.write_html(p.with_suffix(".hist.close.html"))

fig = px.histogram(df, x="SVsize", range_x=(0, 10000), nbins=2000, color="SVtype")
fig.write_html(p.with_suffix(".hist.medium.html"))

fig = px.histogram(df, x="SVsize", nbins=1000, color="SVtype")
fig.write_html(p.with_suffix(".hist.full.html"))
# %%
