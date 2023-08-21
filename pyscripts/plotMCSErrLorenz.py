import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os


# Define a function to extract the number from the file name
def get_number(file_name):
    return int(file_name.split("gauss_")[1].split(".")[0])


def process_err_each_filter(filter_type, folder_path):
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
    truth_df = truth_df.iloc[:, 0:7]

    nan_file_num = 0
    nan_file_names = []

    # Create the plot
    fig, ax = plt.subplots(5, 1, figsize=(5, 4))
    # Loop through each file
    for file_name in trial_file_names:
        # Create an empty dataframe to store the data
        df = pd.DataFrame()
        df["tSec"] = truth_df["tSec"]

        file_path = os.path.join(folder_path, file_name)

        # Read the trial CSV file into a pandas dataframe
        trial_df = pd.read_csv(file_path)
        trial_df = trial_df.iloc[:, 0:6]

        # Check which elements are NaN using isna()
        nan_mask = trial_df.isna()
        # Count the number of NaN values in the entire dataframe
        total_num_nan = nan_mask.sum().sum()
        print(total_num_nan)

        if total_num_nan > 0:
            # Count the files that contain NaN
            print("running to here")
            nan_file_num += 1
            nan_file_names.append(file_path)
            print(file_path)
        else:
            # Add the errors to the empty dataframe
            df["x1"] = trial_df["EST X1"] - truth_df["x1"]
            df["x2"] = trial_df["EST X2"] - truth_df["x2"]
            df["x3"] = trial_df["EST X3"] - truth_df["x3"]
            df["x4"] = trial_df["EST X4"] - truth_df["x4"]
            df["x5"] = trial_df["EST X5"] - truth_df["x5"]

            if (np.abs(df["x1"]) > 1e7).any():
                continue
            else:
                segments = file_name.split("_")
                err_file_name = folder_path + "err_" + file_name + ".csv"
                # save dataframe to csv
                df.to_csv(err_file_name, index=False)

                label = segments[1]
                ax[0].plot(df["tSec"], df["x1"], linewidth=1, label=label)
                ax[0].set_ylabel("x1")
                ax[1].plot(df["tSec"], df["x2"], linewidth=1, label=label)
                ax[1].set_ylabel("x2")
                ax[2].plot(df["tSec"], df["x3"], linewidth=1, label=label)
                ax[2].set_ylabel("x3")
                ax[3].plot(df["tSec"], df["x4"], linewidth=1, label=label)
                ax[3].set_ylabel("x4")
                ax[4].plot(df["tSec"], df["x5"], linewidth=1, label=label)
                ax[4].set_ylabel("x5")

    # # print all NaN files
    # print(nan_file_names)
    fig.suptitle(filter_type + " MCS errors")

    plt.tight_layout()

    start_index = folder_path.index("_") + 1
    end_index = folder_path.rindex("/")
    keyword = folder_path[start_index:end_index]
    # plt.show()
    fig.savefig("plots/" + filter_type + "_" + keyword + "_MCS_err.pdf")


if len(sys.argv) < 2:
    folder_path = "out/out_lorenz/"
else:
    folder_path = sys.argv[1]

filters = ["ukf", "cut4", "cut6"]
# filters = ["cut6"]
# filters = ["ukf", "cut4", "cut6", "house"]
# filters = ["ukf", "cut4", "cut6", "house", "srhouse"]
# filters = ["house", "ukf"]

for filter_type in filters:
    process_err_each_filter(filter_type, folder_path)
