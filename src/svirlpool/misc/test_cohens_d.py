# test cohen's d on different dummy data sets to get an intuition about its behavior
# %%
import numpy as np
import pandas as pd
from tqdm import tqdm


def cohen_d(x: list, y: list) -> float:
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(
        ((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof
    )


# create two different normal distributed sets of numbers around means a and b


def generate_test_data(
    n: int, n_a: int, mean_a_real: float, mean_b_real: float, std_real: float
) -> list[int]:
    n = int(n)
    n_a = int(n_a)
    assert n > 1, "n must be greater than 1"
    assert n_a > 0, "n_a must be greater than 0"
    n_b = n - n_a
    numbers_a = np.random.normal(mean_a_real, std_real, n_a).astype(int)
    numbers_b = np.random.normal(mean_b_real, std_real, n_b).astype(int)
    return sorted([*numbers_a, *numbers_b])


# %%

# def plot_test(n:int, mean_a_real:float, mean_b_real:float, balance:float, std_real:float) -> None:
#     n_a = max(int(round(n * balance*0.5)), 1)
#     numbers_a, numbers_b, d = test(n, n_a, mean_a_real, mean_b_real, std_real)
#     combined_numbers = [(number,'skyblue') for number in numbers_a]
#     combined_numbers += [(number,'darkorange') for number in numbers_b]
#     # make a scatter plot of the sorted numbers
#     combined_numbers = sorted(combined_numbers, key=lambda x: x[0])
#     plt.figure(figsize=(6,3))
#     plt.scatter(range(len(combined_numbers)), [x[0] for x in combined_numbers], c=[x[1] for x in combined_numbers])
#     plt.title(f"cohen's d: {d}")
#     plt.show()

# %%

params_n = [30]
params_std = [int(round(5 ** (i / 5))) for i in range(5, 16)]
# params_na = range(3, 28, 3) # should follow a normal distribution
params_na_normal = list(
    map(lambda x: max(x, 0), sorted(np.random.normal(15, 5, 20).astype(int)))
)
params_means = [100 + int(round(50 ** (i / 10))) for i in range(6, 15)]
params_hom = np.array(
    np.meshgrid(params_n, params_na_normal, params_means, params_means, params_std)
).T.reshape(-1, 5)

params_means_a = [int(mean / 2) for mean in params_means]
params_means_b = [int(1.5 * mean) for mean in params_means]
params_het = np.array(
    np.meshgrid(params_n, params_na_normal, params_means_a, params_means_b, params_std)
).T.reshape(-1, 5)

params = np.concatenate([params_hom, params_het])

print(f"{len(params)} parameter sets to test")
# %%
# concatenate the hom and het parameter sets

# create test dataframe with features and classification (hom, het)
generated_data = []
results = np.zeros(len(params))
# then compute the results
for i, (n, na, mean_a, mean_b, std) in tqdm(enumerate(params)):
    numbers = generate_test_data(n, na, mean_a, mean_b, std)
    generated_data.append(numbers)
    d = cohen_d(numbers[:na], numbers[na:])
    results[i] = d

# %%
df = pd.DataFrame(params, columns=["n", "na", "mean_a", "mean_b", "std"])
df["mean_total"] = (df["mean_a"] + df["mean_b"]) / 2
df["cohens_d"] = results
df["class"] = ["hom"] * len(params_hom) + ["het"] * len(params_het)
df["index"] = range(len(generated_data))

# %%
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

X = df[["n", "cohens_d", "std", "mean_total"]]
Y = df["class"]
X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.20, random_state=42
)
clf = RandomForestClassifier(random_state=42)
clf.fit(X_train, Y_train)
Y_pred = clf.predict(X_test)
# score
print(accuracy_score(Y_test, Y_pred))

# subset df to the test sub set
df_test = df.loc[X_test.index]
# maka column pred_result and classes False or True
df_test["pred_result"] = [
    "True" if Y_pred[i] == Y_test.iloc[i] else "False" for i in range(len(Y_test))
]
# %%
# visualize with a scatter plot
# a dot is a test case with axis x = mean_a, y=mean_b, color = cohens_d, symbol=
