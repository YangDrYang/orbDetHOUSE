import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import processing
import seaborn as sns
import numpy as np
import sys

# filters = ["ukf", "cut4", "cut6", "house", "srhouse"]
# filters = ["ukf", "cut4", "cut6", "house"]
# filters = ["house", "ukf", "cut4"]
# filters = ["house", "cut4", "cut6"]
filters = ["house", "srhouse", "ukf", "srukf"]
# filters = ["srhouse"]

if len(sys.argv) < 2:
    folder_path = "out/out_dense/"
else:
    folder_path = sys.argv[1]
start_index = folder_path.index("_") + 1
end_index = folder_path.rindex("/")
keyword = folder_path[start_index:end_index]

# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(ncols=2, figsize=(10, 5))

for filter_type in filters:
    processing.process_rmse_each_filter(filter_type, folder_path)

    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "trajectory_error_" + filter_type + ".csv")

    # Plot the data on the axes
    ax[0].semilogy(df["tSec"], df["pos_err_rms"], label=filter_type)

    # Plot the data on the axes
    ax[1].semilogy(df["tSec"], df["vel_err_rms"], label=filter_type)

# Add labels and title
ax[0].set_xlabel("time elapsed (s)")
ax[0].set_ylabel("position errors (m)")
# ax_pos.set_ylim([-1000, 5000])
ax[0].set_title("3D position rmse")
ax[0].legend()

# Add labels and title
ax[1].set_xlabel("time elapsed (s)")
ax[1].set_ylabel("velocity errors (m/s)")
# ax[1].set_ylim([-0.5, 1])
ax[1].set_title("3D velocity rmse")
ax[1].legend()

plt.tight_layout()
fig.savefig("plots/all_rmse_" + keyword + ".pdf")

# Create a new figure and axes for violin plot
fig, ax = plt.subplots(ncols=2, figsize=(10, 5))
pos_rmse = []
vel_rmse = []
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "trajectory_error_" + filter_type + ".csv")

    pos_rmse.append(df["pos_err_rms"])
    vel_rmse.append(df["vel_err_rms"])


# Create the violin plot for pos errors
sns.violinplot(
    data=pos_rmse, inner="box", linewidth=1, ax=ax[0], scale="count", scale_hue=False
)

# Add x and y axis labels and a title
ax[0].set_xlabel("filters")
ax[0].set_ylabel("3d pos rmse")
ax[0].yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: "$10^{{{}}}$".format(int(x)))
)
ax[0].set_yscale("log")
ax[0].set_title("violin plot of 3D pos rmse")
ax[0].set_xticklabels(filters, fontsize=12)

# Add annotations for median values
medians = [round(np.median(x), 2) for x in pos_rmse]
for i in range(len(medians)):
    ax[0].annotate(
        "Median: {}".format(medians[i]),
        xy=(i, -4.5),
        ha="center",
        va="center",
        fontsize=12,
    )

# Create the violin plot for pos errors
sns.violinplot(
    data=vel_rmse, inner="box", linewidth=1, ax=ax[1], scale="count", scale_hue=False
)

ax[1].set_xlabel("filters")
ax[1].set_ylabel("3d vel rmse")
ax[1].yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: "$10^{{{}}}$".format(int(x)))
)
ax[1].set_yscale("log")
ax[1].set_title("violin plot of 3D vel rmse")
ax[1].set_xticklabels(filters, fontsize=12)

# Add annotations for median values
medians = [round(np.median(x), 2) for x in vel_rmse]
for i in range(len(medians)):
    ax[1].annotate(
        "Median: {}".format(medians[i]),
        xy=(i, -4.5),
        ha="center",
        va="center",
        fontsize=12,
    )

plt.tight_layout()
fig.savefig("plots/all_rmse_violin_" + keyword + ".pdf")
