import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Define a function to extract the number from the file name
def get_number(file_name):
    return int(file_name.split("_")[1].split(".")[0])


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

        if total_num_nan > 2000:
            # Count the files that contain NaN
            nan_file_names.append(file_path)
            # print(file_path)
        else:
            if (np.abs(trial_df["EST X1"] - truth_df["x"]) > 1e7).any():
                large_value_file_names.append(file_path)
            else:
                # Add the errors to the empty dataframe
                df["pos_err_x"] = (
                    df.get("pos_err_x", 0) + trial_df["EST X1"] - truth_df["x"]
                )
                df["pos_err_y"] = (
                    df.get("pos_err_y", 0) + trial_df["EST X2"] - truth_df["y"]
                )
                df["pos_err_z"] = (
                    df.get("pos_err_z", 0) + trial_df["EST X3"] - truth_df["z"]
                )
                df["vel_err_x"] = (
                    df.get("vel_err_x", 0) + trial_df["EST X4"] - truth_df["vx"]
                )
                df["vel_err_y"] = (
                    df.get("vel_err_y", 0) + trial_df["EST X5"] - truth_df["vy"]
                )
                df["vel_err_z"] = (
                    df.get("vel_err_z", 0) + trial_df["EST X6"] - truth_df["vz"]
                )

    # print all NaN files
    print(nan_file_names)
    # print all large value files
    print(large_value_file_names)
    # average of all trials
    df = df.div(
        len(trial_file_names) - len(nan_file_names) - len(large_value_file_names)
    )

    df["tSec"] = truth_df["tSec"]
    # df["time_lapse"] = truth_df["t"]
    # define the new order of the columns
    new_order = [
        "tSec",
        "pos_err_x",
        "pos_err_y",
        "pos_err_z",
        "vel_err_x",
        "vel_err_y",
        "vel_err_z",
    ]
    df = df.reindex(columns=new_order)

    # calculate the rms of position and velocity for each epoch
    # select the position components to include in the calculation
    cols_pos_err = ["pos_err_x", "pos_err_y", "pos_err_z"]
    # select the velocity components to include in the calculation
    cols_vel_err = ["vel_err_x", "vel_err_y", "vel_err_z"]

    # calculate the root mean square of the selected columns
    rms_pos_err = np.sqrt(np.mean(np.square(df[cols_pos_err]), axis=1))
    rms_vel_err = np.sqrt(np.mean(np.square(df[cols_vel_err]), axis=1))

    # add the root mean square as a new column to the original DataFrame
    df["pos_err_rms"] = rms_pos_err
    df["vel_err_rms"] = rms_vel_err

    # # Display the final dataframe with the added columns
    print(df)
    # save the errors
    df.to_csv("plots/" + filter_type + "_trajectory_error.csv")


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
    fig, ax = plt.subplots(3, 2)
    # Loop through each file
    for file_name in trial_file_names:
        # Create an empty dataframe to store the data
        df = pd.DataFrame()
        df["tSec"] = truth_df["tSec"]

        file_path = os.path.join(folder_path, file_name)

        # Read the trial CSV file into a pandas dataframe
        trial_df = pd.read_csv(file_path)
        trial_df = trial_df.iloc[:, 0:7]

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
            # print(file_path)
        else:
            # Add the errors to the empty dataframe
            df["pos_err_x"] = trial_df["EST X1"] - truth_df["x"]
            df["pos_err_y"] = trial_df["EST X2"] - truth_df["y"]
            df["pos_err_z"] = trial_df["EST X3"] - truth_df["z"]
            df["vel_err_x"] = trial_df["EST X4"] - truth_df["vx"]
            df["vel_err_y"] = trial_df["EST X5"] - truth_df["vy"]
            df["vel_err_z"] = trial_df["EST X6"] - truth_df["vz"]

            if (np.abs(df["pos_err_x"]) > 1e7).any():
                continue
            else:
                segments = file_name.split("_")
                err_file_name = folder_path + "err_" + file_name + ".csv"
                # save dataframe to csv
                df.to_csv(err_file_name, index=False)

                label = segments[1]
                ax[0, 0].plot(df["tSec"], df["pos_err_x"], label=label)
                ax[0, 0].set_ylim(-60000, 20000)
                ax[1, 0].plot(df["tSec"], df["pos_err_y"], label=label)
                ax[1, 0].set_ylim(-60000, 20000)
                ax[2, 0].plot(df["tSec"], df["pos_err_z"], label=label)
                ax[2, 0].set_ylim(-60000, 20000)
                ax[0, 1].plot(df["tSec"], df["vel_err_x"], label=label)
                ax[1, 1].plot(df["tSec"], df["vel_err_y"], label=label)
                ax[2, 1].plot(df["tSec"], df["vel_err_z"], label=label)

    # # print all NaN files
    # print(nan_file_names)
    fig.suptitle(filter_type + " MCS errors")
    # plt.show()
    fig.savefig("plots/" + filter_type + "_MCS_err.pdf")


# # filter_type = "ukf"
# filter_type = "house"
# process_each_filter(filter_type)
