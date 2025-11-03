# %%
import typing
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from tqdm import tqdm

# %%


# function to find if a tuple (chr, start, end) is within a range of 30 of and interval (chr, start, end) in df_repeats
# and returns the chr, start and end of the repeat
def is_in_interval(
    chr: str, start: int, end: int, df_repeats: pd.DataFrame
) -> typing.Tuple[str, int, int]:
    # get the interval from the df_repeats
    interval = df_repeats[df_repeats["chr"] == chr]
    # x = start
    # y = end
    # a = interval['start']
    # b = interval['end']
    interval = interval[
        (end <= interval["end"]) & (end >= interval["start"])
        | (interval["start"] >= start) & (interval["start"] <= end)
        | (start >= interval["start"]) & (start <= interval["end"])
        | (interval["end"] <= end) & (interval["end"] >= start)
    ]
    # if the interval is empty return None
    if interval.empty:
        return None
    # else return the chr, start and end of the interval
    else:
        return interval.iloc[0, 0], interval.iloc[0, 1], interval.iloc[0, 2]


# function to get a np.array from a given memmap in the interval (chr, start, end) which is offset by the chromosome
def get_interval_from_memmap(
    chr: str, start: int, end: int, memmap: np.memmap, dict_memmap_index: dict
) -> np.ndarray:
    # the real index in the memmap is offset by dict_memmap_index[chr]
    offset = dict_memmap_index[chr]
    # get the interval from the memmap
    return np.array(memmap[offset + start : offset + end], dtype=np.float16)


# %%
# =============================================================================
#  LOAD DATA
# =============================================================================
# load signal data
path_fai = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/hs37d5.fa.fai"
)
df_fai = pd.read_csv(path_fai, sep="\t", header=None)
# crate dictionary with key=i and value=df_fai.iloc[i,0]
dict_fai = dict(zip(range(len(df_fai)), df_fai.iloc[:, 0]))

size_memmap = df_fai.iloc[:, 1].sum()
df_memmap_index = df_fai.iloc[:, 1].cumsum()
# add 0 to the beginning of the index and reindex the df
df_memmap_index = pd.concat([pd.Series([0]), df_memmap_index])
df_memmap_index = df_memmap_index.reset_index()

# make dictionary from df_memmap_indexm with key=df_fai.iloc[i,0] and value df_memmap_index.iloc[i]
dict_memmap_index = dict(zip(df_fai.iloc[:, 0], df_memmap_index.iloc[:, 1]))

# open memmap complexity
path_memmap_complexity = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/hs37d5.complexity.npy"
)
memmap_complexity = np.memmap(
    path_memmap_complexity, dtype=np.float16, mode="r", shape=size_memmap
)

# open memmap GC
path_memmap_GC = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/hs37d5.GC.npy"
)
memmap_GC = np.memmap(path_memmap_GC, dtype=np.float16, mode="r", shape=size_memmap)

# load repeats
path_repeats = Path(
    "/home/memsonmi/development/LRSV-detection/tools/pbsv/annotations/human_hs37d5.trf.bed"
)
# read repeats as tsv df
df_repeats = pd.read_csv(path_repeats, sep="\t", header=None)
df_repeats.columns = ["chr", "start", "end"]
df_repeats = df_repeats.astype({"chr": str, "start": int, "end": int})

# load truth be file with indels
path_truth = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/truth.bed"
)
df_truth = pd.read_csv(path_truth, sep="\t", header=None)
df_truth.columns = ["chr", "start", "end", "c0"]
# split c0 column into svtype and svsize. c0 is fomratted as "DEL:1000"
df_truth[["svtype", "svsize"]] = df_truth["c0"].str.split(":", expand=True)
# drop c0 column
df_truth = df_truth.drop(columns=["c0"])
# set correct data types
df_truth = df_truth.astype(
    {"chr": str, "start": int, "end": int, "svtype": str, "svsize": int}
)

path_signal = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/ont-r10.32x/minimap2.ont-r10.20x.signal"
)
# load tsv file with signal data
df_signal = pd.read_csv(path_signal, sep="\t", header=None)
df_signal.columns = ["chr", "start", "end", "svtype", "svsize", "readID"]
df_signal["chr"] = df_signal["chr"].map(dict_fai)
# %%
# =============================================================================
#  PROCESS DATA
# =============================================================================

# use bedtools to get the signal data for each indel in the truth set
# truth bed file already exists
# create a bed file from the signal set
path_bed_signal = Path(
    "/home/memsonmi/development/LRSV-detection/development/test/ont-r10.32x/signal.bed"
)
df_signal.to_csv(path_bed_signal, sep="\t", header=False, index=False)

# %%
L = []
# for each signal in truth, find nearby signals
for i in tqdm(range(0, df_truth.shape[0], 18)):
    row = df_truth.iloc[i, :]
    chr, start, end, svtype, svsize = row
    svtype = svtype.upper()
    r_far = 200
    repeat_interval = is_in_interval(chr, start - 10, end + 10, df_repeats)
    window = (start - r_far, end + r_far)
    in_repeat = repeat_interval is not None
    if in_repeat:
        # update the window
        window = (repeat_interval[1], repeat_interval[2])
    # get all signals in the interval (chr, start-r, end+r) with matching svtype
    dft = df_signal[df_signal["chr"] == chr]
    dft = dft[dft["svtype"] == svtype]
    dft = dft[(dft["start"] >= window[0]) & (dft["end"] <= window[1])]
    # dft svsize must be greater than svsize*0.9 and less than svsize*1.1
    mask = (dft["svsize"] >= svsize * 0.9) & (dft["svsize"] <= svsize * 1.1)
    dft_signal = dft[mask]
    dft_noise = dft[~mask]
    noise_signals = dft_noise.shape[0]
    noise_bp = dft_noise["svsize"].sum()
    # compute std of start, end and svsize
    std_start = dft_signal["start"].std()
    std_end = dft_signal["end"].std()
    max_std_position = max(std_start, std_end)
    std_svsize = dft_signal["svsize"].std()
    number_of_signals = dft_signal.shape[0]
    # compute GC and complexity
    # get interval from memmap complexity
    r = 30
    if svtype == "DEL":
        # get intervals from memmap_complexity both for start and end with radius r
        complexity_start_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, start - r, start + r, memmap_complexity, dict_memmap_index
            ),
            0.25,
        )
        complexity_end_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, end - r, end + r, memmap_complexity, dict_memmap_index
            ),
            0.25,
        )
        # calc mean of complexity_start_Q1 and complexity_end_Q1
        complexity_Q1: float = (complexity_start_Q1 + complexity_end_Q1) / 2
        # get intervals from memmap_GC both for start and end with radius r
        GC_content_start_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, start - r, start + r, memmap_GC, dict_memmap_index
            ),
            0.25,
        )
        GC_content_end_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, end - r, end + r, memmap_GC, dict_memmap_index
            ),
            0.25,
        )
        # calc mean of GC_content_start_Q1 and GC_content_end_Q1
        GC_content_Q1: float = (GC_content_start_Q1 + GC_content_end_Q1) / 2
        # is signal close to or in a repeat region?*
    if svtype == "INS":
        complexity_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, start - r, end + r, memmap_complexity, dict_memmap_index
            ),
            0.25,
        )
        GC_content_Q1 = np.quantile(
            get_interval_from_memmap(
                chr, start - r, end + r, memmap_GC, dict_memmap_index
            ),
            0.25,
        )
    # add data to L
    L.append(
        [
            chr,
            start,
            end,
            svtype,
            svsize,
            number_of_signals,
            max_std_position,
            std_svsize,
            complexity_Q1,
            GC_content_Q1,
            in_repeat,
            noise_signals,
            noise_bp,
        ]
    )

# %%
df_analysis = pd.DataFrame(
    L,
    columns=[
        "chr",
        "start",
        "end",
        "svtype",
        "svsize",
        "number_of_signals",
        "max_std_position",
        "std_svsize",
        "complexity_Q1",
        "GC_content_Q1",
        "in_repeat",
        "noise_signals",
        "noise_bp",
    ],
)
# %%
df_analysis.to_csv(
    "/home/memsonmi/development/LRSV-detection/development/test/ont-r10.32x/analysis.tsv",
    sep="\t",
    index=False,
)
# %%
df_analysis = pd.read_csv(
    "/home/memsonmi/development/LRSV-detection/development/test/ont-r10.32x/analysis.tsv",
    sep="\t",
)
# %%
# investigate corelation between complexity, GC content, and in_repeat with max_std_position
df_analysis.loc[
    :, ["complexity_Q1", "GC_content_Q1", "in_repeat", "max_std_position"]
].corr()
# %%
px.scatter(
    df_analysis,
    x="complexity_Q1",
    y="max_std_position",
    facet_col="in_repeat",
    color="in_repeat",
    opacity=0.6,
    log_y=True,
)
# %%
px.scatter(
    df_analysis,
    x="GC_content_Q1",
    y="max_std_position",
    facet_col="in_repeat",
    color="in_repeat",
    opacity=0.6,
    log_y=True,
)
# %%
px.scatter(
    df_analysis,
    x="complexity_Q1",
    y="noise_signals",
    facet_col="in_repeat",
    color="in_repeat",
    opacity=0.6,
    log_y=True,
)
# %%
px.scatter(
    df_analysis,
    x="GC_content_Q1",
    y="noise_signals",
    facet_col="in_repeat",
    color="in_repeat",
    opacity=0.6,
    log_y=True,
)
# %%
df_analysis.corr()

# %%
pca = PCA(n_components=5)
pca.fit(
    df_analysis.loc[
        :,
        [
            "max_std_position",
            "complexity_Q1",
            "GC_content_Q1",
            "noise_signals",
            "in_repeat",
        ],
    ].dropna()
)
print(pca.explained_variance_ratio_)
# plot contributing varibales of pca in a bar plot
px.bar(pca.explained_variance_ratio_)
# %%
# random forest for feature importance on a regression problem
from matplotlib import pyplot
from sklearn.ensemble import RandomForestRegressor

# define dataset
df = df_analysis.loc[
    :,
    [
        "max_std_position",
        "complexity_Q1",
        "GC_content_Q1",
        "noise_signals",
        "in_repeat",
    ],
].dropna()
# split df into train and test data
df_train, df_test = train_test_split(df, test_size=0.2)
X = df_train.loc[
    :, ["complexity_Q1", "GC_content_Q1", "noise_signals", "in_repeat"]
].to_numpy()
y = df_train.loc[:, ["max_std_position"]].to_numpy().ravel()
# define the model
model = RandomForestRegressor()
# fit the model
model.fit(X, y)
# get importance
importance = model.feature_importances_
# summarize feature importance
for i, v in enumerate(importance):
    print(
        "Feature: %s, Score: %.5f"
        % (["complexity_Q1", "GC_content_Q1", "noise_signals", "in_repeat"][i], v)
    )
# plot feature importance
pyplot.bar([x for x in range(len(importance))], importance)
pyplot.show()
# %%
# validate model
X_test = df_test.loc[
    :, ["complexity_Q1", "GC_content_Q1", "noise_signals", "in_repeat"]
].to_numpy()
y_test = df_test.loc[:, ["max_std_position"]].to_numpy().ravel()
y_pred = model.predict(X_test)
# %%
