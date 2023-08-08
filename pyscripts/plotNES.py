import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import processing
import seaborn as sns
import numpy as np
import scipy.stats as stats

# reference: https://kalman-filter.com/normalized-estimation-error-squared/

filters = ["ukf", "cut4", "cut6", "house", "srhouse"]
# filters = ["house", "ukf"]
# filters = ["cut4"]

folder_path = "out/out_dense/"
# folder_path = "out/out_sparse/"
start_index = folder_path.index("_") + 1
end_index = folder_path.rindex("/")
keyword = folder_path[start_index:end_index]

# flag to determine to process the NES or use existing files for plots
flag_proc = 1

# Set the alpha value
alpha = 0.05
# Calculate the degrees of freedom
dof = 6 - 1
# Calculate the chi-square value at the alpha/2 probability level
chi2_lower = stats.chi2.ppf(alpha / 2, dof)
# Calculate the chi-square value at the 1-alpha/2 probability level
chi2_upper = stats.chi2.ppf(1 - alpha / 2, dof)

# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(ncols=1, figsize=(5, 4))

for filter_type in filters:
    if flag_proc:
        processing.process_nes_each_filter(filter_type, folder_path)

    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "nes_" + filter_type + ".csv")

    # Plot the data on the axes
    ax.scatter(df["tSec"], df["NES"], s=3, label=filter_type)

    if filter_type == filters[-1]:
        # Add a horizontal line for the upper bound
        ax.axhline(y=chi2_upper, color="r", linestyle="--", label="upper bound")
        # Add a horizontal line for the lower bound
        ax.axhline(y=chi2_lower, color="k", linestyle="--", label="lower bound")
    else:
        # Add a horizontal line for the upper bound
        ax.axhline(y=chi2_upper, color="r", linestyle="--")
        # Add a horizontal line for the lower bound
        ax.axhline(y=chi2_lower, color="k", linestyle="--")

# Add labels and title
ax.set_xlabel("time elapsed (s)")
ax.set_ylabel("normalised error square (m)")
ax.set_ylim([-1, 15])
# ax.set_title("nes")
ax.legend()


plt.tight_layout()
fig.savefig("plots/all_nes_." + keyword + ".pdf")
