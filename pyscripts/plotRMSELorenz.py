import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os


# Define a function to extract the number from the file name
def get_number(file_name):
    return int(file_name.split("gauss_")[1].split(".")[0])


def process_rmse_each_filter(filter_type, folder_path):
    # folder_path = "out_/"  # replace with the path to your folder
    # Get a sorted list of file names in the folder that start with filter_type

    # Sort the file names based on the extracted number
    trial_file_names = sorted(
        [
            filename
            for filename in os.listdir(folder_path)
            if filename.startswith(filter_type) and filename.endswith(".csv")
        ],
        key=get_number,
    )

    # print(trial_file_names)

    # Read the truth CSV file into a pandas dataframe
    truth_df = pd.read_csv(folder_path + "trajectory_truth.csv")
    # print(truth_df.head())
    truth_df = truth_df.iloc[:, 0:7]

    # Create an empty dataframe to store the data
    df = pd.DataFrame()
    # df = pd.DataFrame({"t": truth_df.loc[:, "t"]})
    # df = df.rename(columns={"t": "time_lapse"})
    # print(trial_file_names[:20])
    nan_file_names = []
    large_value_file_names = []
    # Loop through each file
    for file_name in trial_file_names:
        file_path = os.path.join(folder_path, file_name)

        # Read the trial CSV file into a pandas dataframe
        trial_df = pd.read_csv(file_path)
        trial_df = trial_df.iloc[:, 0:7]
        # print(trial_df.head())

        # Check which elements are NaN using isna()
        nan_mask = trial_df.isna()
        # Count the number of NaN values in the entire dataframe
        total_num_nan = nan_mask.sum().sum()
        # print(total_num_nan)

        if total_num_nan > 10:
            # Count the files that contain NaN
            nan_file_names.append(file_path)
            # print(file_path)
        else:
            if (np.abs(trial_df["EST X1"] - truth_df["x1"]) > 1e7).any():
                large_value_file_names.append(file_path)
            else:
                # Add the errors to the empty dataframe
                df["x1"] = df.get("x1", 0) + trial_df["EST X1"] - truth_df["x1"]
                df["x2"] = df.get("x2", 0) + trial_df["EST X2"] - truth_df["x2"]
                df["x3"] = df.get("x3", 0) + trial_df["EST X3"] - truth_df["x3"]
                df["x4"] = df.get("x4", 0) + trial_df["EST X4"] - truth_df["x4"]
                df["x5"] = df.get("x5", 0) + trial_df["EST X5"] - truth_df["x5"]

    # print all NaN files
    print("nan file names by " + filter_type + ":   ")
    print(nan_file_names)
    # print all large value files
    print("large value file names by " + filter_type + ":   ")
    print(large_value_file_names)

    # average of all trials
    df = df.div(
        len(trial_file_names) - len(nan_file_names) - len(large_value_file_names)
    )
    # df = df.div(len(trial_file_names) - len(nan_file_names))

    df["tSec"] = truth_df["tSec"]
    # df["time_lapse"] = truth_df["t"]
    # define the new order of the columns
    new_order = [
        "tSec",
        "x1",
        "x2",
        "x3",
        "x4",
        "x5",
    ]
    df = df.reindex(columns=new_order)

    # calculate the rms of position and velocity for each epoch
    # select the position components to include in the calculation
    cols_state_err = ["x1", "x2", "x3", "x4", "x5"]

    # calculate the root mean square of the selected columns
    rms_state_err = np.sqrt(np.mean(np.square(df[cols_state_err]), axis=1))

    # add the root mean square as a new column to the original DataFrame
    df["state_err_rms"] = rms_state_err

    # # Display the final dataframe with the added columns
    print(df)
    # save the errors
    df.to_csv(folder_path + "trajectory_error_" + filter_type + ".csv", index=False)


# filters = ["ukf", "cut4", "cut6", "house", "srhouse"]
# filters = ["ukf", "cut4", "cut6", "house"]
# filters = ["house", "ukf", "cut4"]
# filters = ["house", "cut4", "cut6"]
# filters = ["house", "ukf"]
# filters = ["srhouse"]
filters = ["ukf", "cut4", "cut6"]

if len(sys.argv) < 2:
    folder_path = "out/out_lorenz/"
else:
    folder_path = sys.argv[1]
start_index = folder_path.index("_") + 1
end_index = folder_path.rindex("/")
keyword = folder_path[start_index:end_index]

# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(figsize=(10, 5))


for filter_type in filters:
    process_rmse_each_filter(filter_type, folder_path)

    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "trajectory_error_" + filter_type + ".csv")

    print(df.head())

    # Plot the data on the axes
    ax.semilogy(df["tSec"], df["state_err_rms"], label=filter_type)

# Add labels and title
ax.set_xlabel("time elapsed (s)")
ax.set_ylabel("state errors (m)")
# ax_pos.set_ylim([-1000, 5000])
ax.set_title("state rmse")
ax.legend()

plt.tight_layout()
fig.savefig("plots/all_rmse_" + keyword + ".pdf")

# Create a new figure and axes for violin plot
fig, ax = plt.subplots(figsize=(10, 5))
state_rmse = []
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "trajectory_error_" + filter_type + ".csv")

    state_rmse.append(df["state_err_rms"])


# Create the violin plot for pos errors
sns.violinplot(
    data=state_rmse, inner="box", linewidth=1, ax=ax, scale="count", scale_hue=False
)

# Add x and y axis labels and a title
ax.set_xlabel("filters")
ax.set_ylabel("state rmse")
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: "$10^{{{}}}$".format(int(x)))
)
ax.set_yscale("log")
ax.set_title("violin plot of state rmse")
ax.set_xticklabels(filters, fontsize=12)

# Add annotations for median values
medians = [round(np.median(x), 2) for x in state_rmse]
for i in range(len(medians)):
    ax.annotate(
        "Median: {}".format(medians[i]),
        xy=(i, -4.5),
        ha="center",
        va="center",
        fontsize=12,
    )
