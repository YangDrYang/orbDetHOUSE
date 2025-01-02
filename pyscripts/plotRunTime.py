import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import processing
import seaborn as sns
import numpy as np

filters = ["house", "ukf", "cut4", "cut6"]

folder_path = "out_sparse/"
# folder_path = "out/"
# folder_path = "out_dense/"

fig, ax = plt.subplots(ncols=1, figsize=(5, 4))
df = pd.DataFrame(columns=filters)
print(df.head())
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df[filter_type] = pd.read_csv(
        folder_path + "run_times_" + filter_type + ".csv",
        usecols=[filter_type],
        squeeze=True,
    )

# Remove zero values from the DataFrame
df = df[df != 0]

# calculate mean of each column
means = df.mean()
print("mean of run times:")
print(means)
# calculate standard deviation of each column
stds = df.std()
print("std of run times:")
print(stds)

ax = df.plot.line(ax=ax)
ax.set_xlabel("Monte Carlo trials")
ax.set_ylabel("run time [s]")
plt.tight_layout()
fig.savefig("plots/run_times.pdf")
