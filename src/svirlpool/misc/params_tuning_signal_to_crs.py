# %%

import subprocess
from math import exp
from shlex import split

import pandas as pd

# import plotly.express as px


# %%


def score_func(results: list):
    return (exp(results[0][0]) * results[0][1] / exp(1.0)) * (
        exp(results[1][0]) * results[1][1] / exp(1.0)
    )


def run_experiments(N: int, tuned_param: str):
    # to each experiment add the parameter dict, the index i, results, and the score
    experiments = []
    for i in range(1, N + 1):
        # 1. set params
        params = {"-s": 0.9, "-k": 50, "-n": 250, "-b": 300}
        params[tuned_param] *= 1.0 + (i * 0.1)
        # make k,n,b values ints
        params["-k"] = int(params["-k"])
        params["-n"] = int(params["-n"])
        params["-b"] = int(params["-b"])
        params_string = " ".join([f"{k} {v}" for k, v in params.items()])
        cmd_execute = f"python3 -m scripts.signal_to_crs \
            -i /home/vinzenz/development/LRSV-detection/development/test/ont-r10.32x.chr16/minimap2.ont-r10.32x.chr16.signal \
            -r /home/vinzenz/development/LRSV-detection/development/test//hs37d5.fa \
            -o /home/vinzenz/development/LRSV-detection/development/test/ont-r10.32x.chr16/minimap2.ont-r10.32x.chr16.crs \
            -R /home/vinzenz/development/LRSV-detection/tools/pbsv/annotations/human_hs37d5.trf.bed \
            {params_string}"
        # 2. execute command
        subprocess.run(split(cmd_execute))

        cmd_bed = f"python3 -m scripts.crs_to_bed \
            -i /home/vinzenz/development/LRSV-detection/development/test/ont-r10.32x.chr16/minimap2.ont-r10.32x.chr16.crs \
            -o /home/vinzenz/development/LRSV-detection/development/test/ont-r10.32x.chr16/QC/pt.{i}.bed"
        subprocess.run(split(cmd_bed))

        cmd_eval = f"python3 -m scripts.QC_crs_truth \
            -i /home/vinzenz/development/LRSV-detection/development/test/ont-r10.32x.chr16/QC/pt.{i}.bed \
            -t /home/vinzenz/development/LRSV-detection/development/test/truth.chr16.bed"

        # execute cmd_eval and read std put to string
        result = subprocess.run(split(cmd_eval), capture_output=True)
        # extract from result the number (float) that comes after each time the string "sensitivity\\t" appears
        sensitivities = [
            float(x.split("sensitivity\t")[1])
            for x in result.stdout.decode("utf-8").split("\n")
            if "sensitivity\t" in x
        ]
        # extract from result the number (float) that comes after each time the string "precision\\t" appears
        precisions = [
            float(x.split("precision\t")[1])
            for x in result.stdout.decode("utf-8").split("\n")
            if "precision\t" in x
        ]
        results = list(zip(sensitivities, precisions))
        score = score_func(results)
        experiments.append([params, i, results, score])
    # evaluate experiments to find the best outcome.
    # best outcome is with sensitivity > 0.98 and maximum precision
    # i,sens,spec,sens,spec
    df = pd.DataFrame([[ex[1], *results[2]] for ex in experiments])
    mask = (df.iloc[1] > 0.98) & (df.iloc[3] > 0.98)
    df = df[mask]
    df["score"] = df.iloc[2] * df.iloc[4]
    df.sort_values(by="score", ascending=False, inplace=True)
    best_index = int(df.iloc[0, 0])
    # return best experiment and all experiments
    return experiments[best_index], experiments


tuned_params = ["-s", "-k", "-n", "-b"]
N = 10
all_experiments = []
counter = 0
for tuned_param in tuned_params:
    counter += 1
    print(counter)
    all_experiments.extend(run_experiments(N, tuned_param))
# %%

all_experiments = [
    [
        {"-s": 0.9900000000000001, "-k": 50, "-n": 250, "-b": 300},
        1,
        [(0.4415, 0.0303), (0.9699, 0.0174)],
        0.00029266187100229367,
    ],
    [
        {"-s": 1.08, "-k": 50, "-n": 250, "-b": 300},
        2,
        [(1.0, 0.1041), (1.0, 0.0019)],
        0.00019778999999999996,
    ],
    [
        {"-s": 1.1700000000000002, "-k": 50, "-n": 250, "-b": 300},
        3,
        [(1.0, 0.1046), (1.0, 0.0019)],
        0.00019873999999999996,
    ],
    [
        {"-s": 1.26, "-k": 50, "-n": 250, "-b": 300},
        4,
        [(1.0, 0.1051), (1.0, 0.0019)],
        0.00019968999999999998,
    ],
    [
        {"-s": 1.35, "-k": 50, "-n": 250, "-b": 300},
        5,
        [(1.0, 0.1055), (1.0, 0.0019)],
        0.00020044999999999994,
    ],
    [
        {"-s": 1.4400000000000002, "-k": 50, "-n": 250, "-b": 300},
        6,
        [(1.0, 0.1057), (1.0, 0.0019)],
        0.00020083,
    ],
    [
        {"-s": 1.5300000000000002, "-k": 50, "-n": 250, "-b": 300},
        7,
        [(1.0, 0.1061), (1.0, 0.0019)],
        0.00020158999999999995,
    ],
    [
        {"-s": 1.62, "-k": 50, "-n": 250, "-b": 300},
        8,
        [(1.0, 0.1065), (1.0, 0.0019)],
        0.00020234999999999996,
    ],
    [
        {"-s": 1.71, "-k": 50, "-n": 250, "-b": 300},
        9,
        [(1.0, 0.1066), (1.0, 0.0019)],
        0.00020253999999999997,
    ],
    [
        {"-s": 1.8, "-k": 50, "-n": 250, "-b": 300},
        10,
        [(1.0, 0.1067), (1.0, 0.0019)],
        0.00020272999999999998,
    ],
    [
        {"-s": 0.9, "-k": 55, "-n": 250, "-b": 300},
        1,
        [(0.4149, 0.0371), (0.9707, 0.0289)],
        0.0005800176677600495,
    ],
    [
        {"-s": 0.9, "-k": 60, "-n": 250, "-b": 300},
        2,
        [(0.4149, 0.037), (0.9707, 0.0289)],
        0.0005784542778199955,
    ],
    [
        {"-s": 0.9, "-k": 65, "-n": 250, "-b": 300},
        3,
        [(0.4149, 0.037), (0.9707, 0.0289)],
        0.0005784542778199955,
    ],
    [
        {"-s": 0.9, "-k": 70, "-n": 250, "-b": 300},
        4,
        [(0.4149, 0.037), (0.9707, 0.0288)],
        0.0005764527059244246,
    ],
    [
        {"-s": 0.9, "-k": 75, "-n": 250, "-b": 300},
        5,
        [(0.4149, 0.037), (0.9707, 0.0288)],
        0.0005764527059244246,
    ],
    [
        {"-s": 0.9, "-k": 80, "-n": 250, "-b": 300},
        6,
        [(0.4149, 0.037), (0.9707, 0.0288)],
        0.0005764527059244246,
    ],
    [
        {"-s": 0.9, "-k": 85, "-n": 250, "-b": 300},
        7,
        [(0.4149, 0.0369), (0.9707, 0.0288)],
        0.0005748947256381424,
    ],
    [
        {"-s": 0.9, "-k": 90, "-n": 250, "-b": 300},
        8,
        [(0.4149, 0.0369), (0.9707, 0.0288)],
        0.0005748947256381424,
    ],
    [
        {"-s": 0.9, "-k": 95, "-n": 250, "-b": 300},
        9,
        [(0.4149, 0.0369), (0.9707, 0.0287)],
        0.0005728985633963433,
    ],
    [
        {"-s": 0.9, "-k": 100, "-n": 250, "-b": 300},
        10,
        [(0.4149, 0.0369), (0.9707, 0.0287)],
        0.0005728985633963433,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 275, "-b": 300},
        1,
        [(0.4149, 0.0371), (0.9707, 0.029)],
        0.0005820246493093923,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 300, "-b": 300},
        2,
        [(0.4149, 0.0371), (0.9707, 0.029)],
        0.0005820246493093923,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 325, "-b": 300},
        3,
        [(0.4149, 0.0372), (0.9707, 0.029)],
        0.0005835934489032181,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 350, "-b": 300},
        4,
        [(0.4149, 0.0372), (0.9707, 0.029)],
        0.0005835934489032181,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 375, "-b": 300},
        5,
        [(0.4149, 0.0372), (0.9707, 0.029)],
        0.0005835934489032181,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 400, "-b": 300},
        6,
        [(0.4149, 0.0372), (0.9707, 0.029)],
        0.0005835934489032181,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 425, "-b": 300},
        7,
        [(0.4202, 0.0377), (0.9727, 0.0291)],
        0.0005978251193553639,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 450, "-b": 300},
        8,
        [(0.4202, 0.0377), (0.9727, 0.0292)],
        0.0005998795012088187,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 475, "-b": 300},
        9,
        [(0.4202, 0.0377), (0.9727, 0.0292)],
        0.0005998795012088187,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 500, "-b": 300},
        10,
        [(0.4202, 0.0378), (0.9727, 0.0292)],
        0.0006014706935197175,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 330},
        1,
        [(0.4149, 0.0378), (0.9707, 0.0273)],
        0.000558243811328501,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 360},
        2,
        [(0.4149, 0.0382), (0.9707, 0.0259)],
        0.0005352203248756637,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 390},
        3,
        [(0.4149, 0.0386), (0.9707, 0.0246)],
        0.0005136790835563031,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 420},
        4,
        [(0.418, 0.0394), (0.971, 0.0238)],
        0.0005090017206421826,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 450},
        5,
        [(0.418, 0.0399), (0.971, 0.0227)],
        0.0004916372994467902,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 480},
        6,
        [(0.418, 0.0403), (0.971, 0.0217)],
        0.0004746908402495364,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 510},
        7,
        [(0.4233, 0.0411), (0.9722, 0.0208)],
        0.0004670615782619177,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 540},
        8,
        [(0.4233, 0.0414), (0.9722, 0.02)],
        0.00045237575659843243,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 570},
        9,
        [(0.4233, 0.0416), (0.9722, 0.0193)],
        0.0004386515065914847,
    ],
    [
        {"-s": 0.9, "-k": 50, "-n": 250, "-b": 600},
        10,
        [(0.4309, 0.0423), (0.9718, 0.0184)],
        0.00042830593886861054,
    ],
]
