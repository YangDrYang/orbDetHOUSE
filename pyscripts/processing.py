import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import datetime as dt
import matplotlib.dates as mdates


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
    df.to_csv(folder_path + "trajectory_error_" + filter_type + ".csv", index=False)


def process_rmse_each_filter_ccdata(
    filter_type, out_folder_path, norad_id, od_ref_data_file, state_type
):
    # folder_path = "out_/"  # replace with the path to your folder
    # Get a sorted list of file names in the folder that start with filter_type

    # Read the truth CSV file into a pandas dataframe
    truth_df = pd.read_csv(od_ref_data_file)
    # print(truth_df.head())
    truth_df = truth_df.iloc[:, 0:7]
    # Read the estimaiton CSV file into a pandas dataframe

    est_file_name = filter_type + "_id_" + str(norad_id) + "_" + state_type + ".csv"
    # Create an empty dataframe to store the data
    est_df = pd.read_csv(out_folder_path + est_file_name)

    df = pd.DataFrame()
    # Add the errors to the empty dataframe
    df["pos_err_x"] = est_df["EST X1"] - truth_df["Interpolated_X"]
    df["pos_err_y"] = est_df["EST X2"] - truth_df["Interpolated_Y"]
    df["pos_err_z"] = est_df["EST X3"] - truth_df["Interpolated_Z"]
    df["vel_err_x"] = est_df["EST X4"] - truth_df["Interpolated_VX"]
    df["vel_err_y"] = est_df["EST X5"] - truth_df["Interpolated_VY"]
    df["vel_err_z"] = est_df["EST X6"] - truth_df["Interpolated_VZ"]

    df["mjd"] = truth_df["MJD"]
    # df["time_lapse"] = truth_df["t"]
    # define the new order of the columns
    new_order = [
        "mjd",
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
    # print(df)
    # save the errors
    df.to_csv(
        out_folder_path
        + filter_type
        + "_err_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".csv",
        index=False,
    )


def process_rmse_HOUSE_ccdata(trial_no, out_folder_path, norad_id, od_ref_data_file):
    # folder_path = "out_/"  # replace with the path to your folder
    # Get a sorted list of file names in the folder that start with filter_type

    # Read the truth CSV file into a pandas dataframe
    truth_df = pd.read_csv(od_ref_data_file)
    # print(truth_df.head())
    truth_df = truth_df.iloc[:, 0:7]
    # Read the estimaiton CSV file into a pandas dataframe

    est_file_name = "house_id_" + str(norad_id) + "_" + str(trial_no) + ".csv"
    # Create an empty dataframe to store the data
    est_df = pd.read_csv(out_folder_path + est_file_name)

    df = pd.DataFrame()
    # Add the errors to the empty dataframe
    df["pos_err_x"] = est_df["EST X1"] - truth_df["Interpolated_X"]
    df["pos_err_y"] = est_df["EST X2"] - truth_df["Interpolated_Y"]
    df["pos_err_z"] = est_df["EST X3"] - truth_df["Interpolated_Z"]
    df["vel_err_x"] = est_df["EST X4"] - truth_df["Interpolated_VX"]
    df["vel_err_y"] = est_df["EST X5"] - truth_df["Interpolated_VY"]
    df["vel_err_z"] = est_df["EST X6"] - truth_df["Interpolated_VZ"]

    df["mjd"] = truth_df["MJD"]
    # df["time_lapse"] = truth_df["t"]
    # define the new order of the columns
    new_order = [
        "mjd",
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
    # print(df)
    # save the errors
    df.to_csv(
        out_folder_path
        + "house_err_id_"
        + str(norad_id)
        + "_"
        + str(trial_no)
        + ".csv",
        index=False,
    )


# Normalised Error Square
def process_nes_each_filter(filter_type, folder_path):
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

    # Initialise an dataframe with two columns
    nes_df = pd.DataFrame()

    nan_file_names = []
    large_value_file_names = []
    # Loop through each file
    for file_name in trial_file_names:
        file_path = os.path.join(folder_path, file_name)

        # Read the trial CSV file into a pandas dataframe
        trial_df = pd.read_csv(file_path)
        # print(trial_df.head())
        trial_state_df = trial_df.iloc[:, 0:7]
        trial_cov_df = trial_df.iloc[:, 7:43]
        # print(trial_cov_df)

        # Check which elements are NaN using isna()
        nan_mask = trial_state_df.isna()
        # Count the number of NaN values in the entire dataframe
        total_num_nan = nan_mask.sum().sum()
        # print(total_num_nan)

        if total_num_nan > 2000:
            # Count the files that contain NaN
            nan_file_names.append(file_path)
            # print(file_path)
        else:
            if (np.abs(trial_state_df["EST X1"] - truth_df["x"]) > 1e7).any():
                large_value_file_names.append(file_path)
            else:
                df = pd.DataFrame(columns=["tSec", "NES"], dtype=float)
                # Iterate over the rows of the dataframes
                for index, row in trial_cov_df.iterrows():
                    # Extract the covariance matrix
                    cov = np.array(row.values.reshape(6, 6))

                    # Extract the state error
                    state_err = np.array(
                        trial_state_df.iloc[index, 1:].values.reshape(6, 1)
                        - truth_df.iloc[index, 1:].values.reshape(6, 1)
                    )

                    # print(state_err)

                    # print(state_err.T @ np.linalg.inv(cov) @ state_err)
                    # Calculate the expression and save to NES column
                    df.loc[index, "NES"] = state_err.T @ np.linalg.inv(cov) @ state_err

                    # print("This is how you pause")

                    # input()
                nes_df["NES"] = nes_df.get("NES", 0) + df["NES"]

    # print all NaN files
    print("nan file names by " + filter_type + ":   ")
    print(nan_file_names)
    # print all large value files
    print("large value file names by " + filter_type + ":   ")
    print(large_value_file_names)

    # average of all trials
    nes_df = nes_df.div(
        len(trial_file_names) - len(nan_file_names) - len(large_value_file_names)
    )
    # nes_df = nes_df.div(len(trial_file_names) - len(nan_file_names))

    nes_df["tSec"] = truth_df["tSec"]

    new_order = [
        "tSec",
        "NES",
    ]
    nes_df = nes_df.reindex(columns=new_order)

    # save the errors
    nes_df.to_csv(folder_path + "nes_" + filter_type + ".csv", index=False)


# Normalised Error Square for cc data
def process_nes_each_filter_ccdata(
    filter_type, out_folder_path, norad_id, od_ref_data_file
):
    # Read the truth CSV file into a pandas dataframe
    truth_df = pd.read_csv(od_ref_data_file)
    truth_df = truth_df.iloc[:, 0:7]

    # print("od reference:    \n", truth_df)

    # Read the estimaiton CSV file into a pandas dataframe

    est_file_name = filter_type + "_id_" + str(norad_id) + ".csv"
    # Create an empty dataframe to store the data
    est_df = pd.read_csv(out_folder_path + est_file_name)
    est_state_df = est_df.iloc[:, 0:7]
    # print("od estimation:   \n", est_state_df)

    est_cov_df = est_df.iloc[:, 7:43]
    # print("od covariance:   \n", est_cov_df)

    # Initialise an dataframe with two columns
    nes_df = pd.DataFrame(columns=["tSec", "NES"], dtype=float)
    nes_df["tSec"] = truth_df["MJD"]
    for index, row in est_cov_df.iterrows():
        # Extract the covariance matrix
        cov = np.array(row.values.reshape(6, 6))

        # Extract the state error
        state_err = np.array(
            est_state_df.iloc[index, 1:].values.reshape(6, 1)
            - truth_df.iloc[index, 1:].values.reshape(6, 1)
        )

        print("od error::   \n", state_err)

        # print(state_err.T @ np.linalg.inv(cov) @ state_err)
        # Calculate the squared NES
        nes_squared = float(state_err.T @ np.linalg.inv(cov) @ state_err)
        # Assign the squared NES to the DataFrame
        nes_df.loc[index, "NES"] = nes_squared

    # save the errors
    nes_df.to_csv(
        out_folder_path + "nes_" + filter_type + "_id_" + str(norad_id) + ".csv",
        index=False,
    )


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
    fig, ax = plt.subplots(3, 2, figsize=(5, 4))
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
            print(file_path)
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
                ax[0, 0].plot(
                    df["tSec"], df["pos_err_x"] / 1000, linewidth=1, label=label
                )
                ax[0, 0].set_ylim(-10, 10)
                ax[0, 0].set_ylabel("x [km]")
                ax[1, 0].plot(
                    df["tSec"], df["pos_err_y"] / 1000, linewidth=1, label=label
                )
                ax[1, 0].set_ylim(-5, 5)
                ax[1, 0].set_ylabel("y [km]")
                ax[2, 0].plot(
                    df["tSec"], df["pos_err_z"] / 1000, linewidth=1, label=label
                )
                ax[2, 0].set_ylim(-10, 10)
                ax[2, 0].set_ylabel("z [km]")
                ax[2, 0].set_xlabel("time elapsed [s]")
                ax[0, 1].plot(df["tSec"], df["vel_err_x"], linewidth=1, label=label)
                ax[0, 1].set_ylim(-5, 5)
                ax[0, 1].set_ylabel("x [m/s]")
                ax[1, 1].plot(df["tSec"], df["vel_err_y"], linewidth=1, label=label)
                ax[1, 1].set_ylim(-3, 3)
                ax[1, 1].set_ylabel("y [m/s]")
                ax[2, 1].plot(df["tSec"], df["vel_err_z"], linewidth=1, label=label)
                ax[2, 1].set_ylim(-5, 5)
                ax[2, 1].set_ylabel("z [m/s]")
                ax[2, 1].set_xlabel("time elapsed [s]")

    # # print all NaN files
    # print(nan_file_names)
    fig.suptitle(filter_type + " MCS errors")

    plt.tight_layout()

    # plt.show()
    fig.savefig("plots/" + filter_type + "_MCS_err.pdf")


def find_closest_row(df, label, identifier):
    closest_row = None
    min_diff = float("inf")

    for index, row in df.iterrows():
        labelled_value = row[label]
        diff = abs(identifier - labelled_value)
        if diff < min_diff:
            closest_row = row
            min_diff = diff

    return closest_row


# Post-residuals for cc data
def process_post_res_each_filter_ccdata(
    od_file, stn_file, meas_file, post_res_file, post_res_plot_file
):
    # Calculate the post residuals

    # print(trial_file_names)

    # Read the real measurement CSV file into a pandas dataframe
    meas_df = pd.read_csv(meas_file)

    # print(meas_df)

    # Read the od result CSV file into a pandas dataframe
    od_df = pd.read_csv(od_file)
    od_df = od_df.iloc[:, 0:7]

    print(od_df)

    # Read the station ECI coordinate CSV file into a pandas dataframe
    stn_df = pd.read_csv(stn_file)

    post_res_df = pd.DataFrame()

    for index, meas_row in meas_df.iterrows():
        # print("meas_row:    ", meas_row)
        # Access the common identifier
        identifier = meas_row["MJD"]
        # print("identifier1: ", identifier)

        stn_row = find_closest_row(stn_df, "MJD", identifier)
        # print("stn_row: ", stn_row)

        post_res_row = pd.DataFrame({"MJD": [identifier]})
        # print(post_res_row)

        # Access the common identifier
        identifier = (meas_row["MJD"] - meas_df.iloc[0, 6]) * 86400

        # print("identifier2: ", identifier)

        # Retrieve the corresponding row in od_df using the common identifier
        od_row = find_closest_row(od_df, "TIME", identifier)

        # Access the values from od_row and meas_row for calculations
        x_sat = od_row["EST X1"]
        y_sat = od_row["EST X2"]
        z_sat = od_row["EST X3"]

        # print("x_sat:")
        # print(x_sat)

        x_stn = stn_row["X_ECI"]
        y_stn = stn_row["Y_ECI"]
        z_stn = stn_row["Z_ECI"]

        # print("x_stn:")
        # print(x_stn)

        # Calculate the range vector between satellite and station
        p = np.vstack([x_sat, y_sat, z_sat]) - np.vstack([x_stn, y_stn, z_stn])
        # print("range vector between satellite and station:    ", p)

        # Calculate residuals of the right ascension angle
        epsilon = 1e-2  # Adjust the epsilon value as needed
        angle_diff = np.arctan2(p[1], p[0]) - meas_row["RA"] / 180 * np.pi
        if angle_diff >= 2 * np.pi - epsilon:
            angle_diff -= 2 * np.pi
        post_res_row["RA"] = angle_diff
        # post_res_row["RA"] = (np.arctan2(p[1], p[0]) - meas_row["RA"] / 180 * np.pi) % (
        #     2 * np.pi
        # )
        # Calculate residuals of the declination angle
        post_res_row["Dec"] = (
            np.arcsin(p[2] / np.linalg.norm(p)) - meas_row["Dec"] / 180 * np.pi
        )

        post_res_df = pd.concat([post_res_df, post_res_row], ignore_index=True)

    print("post_res_df: ", post_res_df)
    # # Generate plots for two angles

    post_res_df.to_csv(post_res_file, index=False)

    # # Convert MJD to datetime objects
    # dates = mdates.num2date(post_res_df["MJD"])

    # # Create subplots with shared x-axis
    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 10))

    # # Plot RA vs MJD
    # ax1.scatter(dates, np.radians(post_res_df["RA"]), color="blue")
    # ax1.set_ylabel("ra res(radians)")
    # # ax1.set_title("ra res vs date")

    # # Plot Dec vs MJD
    # ax2.scatter(dates, np.radians(post_res_df["Dec"]), color="red")
    # ax2.set_xlabel("hour (utc)")
    # ax2.set_ylabel("dec res (radians)")
    # # ax2.set_title("dec vs date")

    # # Format x-axis as hours
    # hours = mdates.HourLocator(interval=2)
    # hour_format = mdates.DateFormatter("%H:%M")
    # ax2.xaxis.set_major_locator(hours)
    # ax2.xaxis.set_major_formatter(hour_format)

    # # Show only the date for the first epoch of each day
    # dates_first_epoch = np.unique([date.date() for date in dates])
    # for date in dates_first_epoch:
    #     ax1.axvline(date, color="gray", linestyle="--", alpha=0.5)
    #     ax2.axvline(date, color="gray", linestyle="--", alpha=0.5)

    # # Show only 1 hour before and 1 hour after the data on x-axis
    # start_time = min(dates) - dt.timedelta(hours=1)
    # end_time = max(dates) + dt.timedelta(hours=1)
    # ax2.set_xlim(start_time, end_time)

    # # Rotate x-axis tick labels and adjust label spacing
    # fig.autofmt_xdate(rotation=45)
    # ax2.tick_params(axis="x", rotation=45, labelsize=10)

    # # Adjust spacing between subplots
    # plt.subplots_adjust(hspace=0.4)

    # plt.tight_layout()
    # # Save the figure
    # plt.savefig(post_res_plot_file)
