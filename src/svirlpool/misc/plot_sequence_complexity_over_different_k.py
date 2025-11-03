# %%
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

# %%
# this bash cmd is used to produce the tsv file used in this script
executable = "/home/vinzenz/development/tutorial/build/fast_sequence_complexity"
f"""
paste \
<(cat hs37d5.only.fa | {executable} -w 15 -k 1 -k 2 -k 3 -k 4 -k 5) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 1) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 2) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 3) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 4) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 5) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 1 -k 2 -k 3) \
<(cat hs37d5.only.fa | {executable} -w 15 -k 2 -k 3 ) > comp.tsv
"""
# %%
path_tsv = Path("/home/vinzenz/development/LRSV-detection/development/test/comp.tsv")
path_fig = Path(
    "/home/vinzenz/development/LRSV-detection/development/test/sequence_complexity.different_k.hs37d5.html"
)
df = pd.read_csv(path_tsv, sep="\t", header=None, dtype=np.float16)

# %%
df.columns = ["12345", *list("12345"), "123", "23"]
# %%
# calculate pairwise differences of all pairs of columns in df
# first, create all 2-combinations of column names
pairs = []
for i in range(len(df.columns)):
    for j in range(i + 1, len(df.columns)):
        pairs.append((df.columns[i], df.columns[j]))

# sample 100000 rows from df and save in a new dataframe df_sample
df_sample = df.sample(100000)

# %%
# then, for each pair, calculate the difference between the two columns
# and add it as a new column to the dataframe
# the new column name is the concatenated pair of column names with a '-' in between
for pair in pairs:
    df_sample[pair[0] + "-" + pair[1]] = (df_sample[pair[0]] - df_sample[pair[1]]).abs()

# %%
# create a plotly express box plot of all columns of df_sample
# set points=False to not show the individual data points
fig = px.box(df_sample, points=False)
# %%
fig.write_html(path_fig)
# %%
