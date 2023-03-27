import os
import pandas as pd
import numpy as np


# Define a function to extract the number from the file name
def get_number(file_name):
    return int(file_name.split("_")[1].split(".")[0])


def process_each_filter(filter_type, folder_path):
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

    # Create an empty dataframe to store the data
    df = pd.DataFrame()
    # df = pd.DataFrame({"t": truth_df.loc[:, "t"]})
    # df = df.rename(columns={"t": "time_lapse"})
    # print(trial_file_names[:20])
    nan_file_num = 0
    nan_file_names = []
    # Loop through each file
    for file_name in trial_file_names:
        file_path = os.path.join(folder_path, file_name)

        # Read the trial CSV file into a pandas dataframe
        trial_df = pd.read_csv(file_path)
        print(trial_df.head())

        # Check which elements are NaN using isna()
        nan_mask = trial_df.isna()
        # Count the number of NaN values in the entire dataframe
        total_num_nan = nan_mask.sum().sum()
        # print(total_num_nan)

        if total_num_nan > 2000:
            # Count the files that contain NaN
            nan_file_num += 1
            nan_file_names.append(file_path)
            # print(file_path)
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
    # average of all trials
    df = df.div(len(trial_file_names) - nan_file_num)

    df["time_lapse"] = truth_df["tSec"]
    # df["time_lapse"] = truth_df["t"]
    # define the new order of the columns
    new_order = [
        "time_lapse",
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


# # filter_type = "ukf"
# filter_type = "house"
# process_each_filter(filter_type)
