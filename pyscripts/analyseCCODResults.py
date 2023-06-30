import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np

filters = ["house", "ukf", "cut4", "cut6"]
# filters = ["ukf"]
# filters = ["house"]
# Directory path
out_folder_path = "out/out_ccdata/"

norad_id = 46984
meas_file = "ccdata/meas_data_id_" + str(norad_id) + ".csv"
stn_file = "ccdata/stn_eci_coordinates.csv"
od_ref_data_file = "refdata/od_ref_id_" + str(norad_id) + ".csv"
post_res_file = out_folder_path + "post_res_id_" + str(norad_id) + ".csv"


meas_df = pd.read_csv(meas_file)
stamp = (meas_df["MJD"] - meas_df["MJD"][0]) * 1440
print(stamp)

# Calculate the time differences between consecutive timestamps
stamp_diff = np.diff(stamp)

# Define a threshold to determine when a new time window starts
threshold = 1000.0  # Adjust the threshold as needed

# Find the indices where the time differences exceed the threshold
window_starts = np.where(stamp_diff > threshold)[0] + 1

# Initialize the time windows list
time_windows = []

# Iterate over the window starts and create the time windows
for i, start_idx in enumerate(window_starts):
    if i == 0:
        prev_idx = 0
    else:
        prev_idx = window_starts[i - 1]

    # Extract the start and end indices for the current window
    start_time = stamp[prev_idx]
    end_time = stamp[start_idx - 1]

    # Append the time window to the list
    time_windows.append((start_time, end_time))

# # Create a new figure and axes for logarithmic plot
# fig, ax = plt.subplots(ncols=2, figsize=(10, 5))
# for filter_type in filters:
#     # post measurement residuals
#     od_file = "out_ccdata/" + filter_type + "_id_" + str(norad_id) + ".csv"
#     # processing.process_post_res_each_filter(od_file, stn_file, meas_file, post_res_file)

#     processing.process_rmse_each_filter_ccdata(
#         filter_type, out_folder_path, norad_id, od_ref_data_file
#     )

#     # Read CSV file into a pandas dataframe
#     err_df = pd.read_csv(
#         out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
#     )

#     stamp = (err_df["mjd"] - err_df["mjd"][0]) * 1440

#     # Plot the data on the axes
#     ax[0].plot(stamp, err_df["pos_err_rms"], label=filter_type)

#     # Plot the data on the axes
#     ax[1].plot(stamp, err_df["vel_err_rms"], label=filter_type)

# # Add labels and title
# ax[0].set_xlabel("time elapsed (s)")
# ax[0].set_ylabel("position errors (m)")
# # ax[0].set_ylim([-10000, 10000])
# ax[0].set_title("3D position rmse")
# ax[0].legend()

# # Add labels and title
# ax[1].set_xlabel("time elapsed (s)")
# ax[1].set_ylabel("velocity errors (m/s)")
# # ax[1].set_ylim([-1, 1])
# ax[1].set_title("3D velocity rmse")
# ax[1].legend()

# plt.tight_layout()
# fig.savefig("plots/all_rmse_id_" + str(norad_id) + ".pdf")

# Create an empty DataFrame to store RMS values
df_pos_rmse = pd.DataFrame(
    columns=["filter_type", "time_window", "x_rmse", "y_rmse", "z_rmse"]
)


# Plot the data with the identified time windows
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))

# Create a color mapping dictionary
color_map = {"house": "orange", "ukf": "red", "cut4": "blue", "cut6": "green"}

for i, filter_type in enumerate(filters):
    processing.process_rmse_each_filter_ccdata(
        filter_type, out_folder_path, norad_id, od_ref_data_file
    )
    # Read the data into a pandas dataframe
    err_df = pd.read_csv(
        out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
    )

    # Get the color for the current filter type from the color mapping dictionary
    color = color_map.get(
        filter_type, "black"
    )  # Use black as the default color if not found in the mapping

    # Plot absolute x errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        x_errors = err_df["pos_err_x"][start_idx]
        y_errors = err_df["pos_err_y"][start_idx]
        z_errors = err_df["pos_err_z"][start_idx]

        # Calculate RMS for each error component within the sub-time window
        x_rmse = np.sqrt(np.mean(x_errors**2))
        y_rmse = np.sqrt(np.mean(y_errors**2))
        z_rmse = np.sqrt(np.mean(z_errors**2))

        # Store the RMS values in the DataFrame
        df_pos_rmse = pd.concat(
            [
                df_pos_rmse,
                pd.DataFrame(
                    {
                        "filter_type": [filter_type],
                        "time_window": [f"{start_time}-{end_time}"],
                        "x_rmse": [x_rmse],
                        "y_rmse": [y_rmse],
                        "z_rmse": [z_rmse],
                    }
                ),
            ],
            ignore_index=True,
        )

        sub_ax = ax[0, j]
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["pos_err_x"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("x absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

    # Plot absolute y errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[1, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["pos_err_y"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("y absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

    # Plot absolute z errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[2, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["pos_err_z"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("z absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

# Set x-axis limits and x-ticks for each subplot to show all three time windows
for i in range(3):
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[i, j]
        sub_ax.set_xlim(start_time, end_time)
        sub_ax.set_xticks([start_time, end_time])

# Set font size for tick labels
for sub_ax in ax.flat:
    sub_ax.tick_params(axis="both", which="both", labelsize=12)
plt.tight_layout()
# plt.show()
fig.savefig("plots/all_pos_abserr_id_" + str(norad_id) + ".pdf")
# Print the DataFrame with RMS values
print(df_pos_rmse)


# Create an empty DataFrame to store RMS values
df_vel_rmse = pd.DataFrame(
    columns=["filter_type", "time_window", "vx_rmse", "vy_rmse", "vz_rmse"]
)

# Plot the data with the identified time windows
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))

for i, filter_type in enumerate(filters):
    # Read the data into a pandas dataframe
    err_df = pd.read_csv(
        out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
    )

    # Get the color for the current filter type from the color mapping dictionary
    color = color_map.get(
        filter_type, "black"
    )  # Use black as the default color if not found in the mapping

    # Plot absolute x errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        vx_errors = err_df["vel_err_x"][start_idx]
        vy_errors = err_df["vel_err_y"][start_idx]
        vz_errors = err_df["vel_err_z"][start_idx]

        # Calculate RMS for each error component within the sub-time window
        vx_rmse = np.sqrt(np.mean(vx_errors**2))
        vy_rmse = np.sqrt(np.mean(vy_errors**2))
        vz_rmse = np.sqrt(np.mean(vz_errors**2))

        # Store the RMS values in the DataFrame
        df_vel_rmse = pd.concat(
            [
                df_vel_rmse,
                pd.DataFrame(
                    {
                        "filter_type": [filter_type],
                        "time_window": [f"{start_time}-{end_time}"],
                        "vx_rmse": [vx_rmse],
                        "vy_rmse": [vy_rmse],
                        "vz_rmse": [vz_rmse],
                    }
                ),
            ],
            ignore_index=True,
        )

        sub_ax = ax[0, j]
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["vel_err_x"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vx absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

    # Plot absolute y errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[1, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["vel_err_y"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vy absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

    # Plot absolute z errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[2, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["vel_err_z"][start_idx]),
            color=color,
            label=filter_type,
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vz absolute error (m)")
        sub_ax.set_yscale("log")
        sub_ax.legend()

# Set x-axis limits and x-ticks for each subplot to show all three time windows
for i in range(3):
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[i, j]
        sub_ax.set_xlim(start_time, end_time)
        sub_ax.set_xticks([start_time, end_time])

# Set font size for tick labels
for sub_ax in ax.flat:
    sub_ax.tick_params(axis="both", which="both", labelsize=12)
plt.tight_layout()
# plt.show()
fig.savefig("plots/all_vel_abserr_id_" + str(norad_id) + ".pdf")
# Print the DataFrame with RMS values
print(df_vel_rmse)

# Define the filter types and time windows
filter_types = ["house", "ukf", "cut4", "cut6"]
time_windows = [
    "0.0-2.1686866565141827",
    "1462.6672116550617-1464.8276216641534",
    "2926.092338307062-2927.989939986728",
]

# Set the figure size and create subplots
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(11, 6))

# Set the font size
plt.rcParams.update({"font.size": 12})

# Calculate the maximum RMSE value across all filter types and time windows
max_pos_rmse = df_pos_rmse[["x_rmse", "y_rmse", "z_rmse"]].max().max() * 1.05

# Loop over the filter types
for i, filter_type in enumerate(filter_types):
    # Calculate the subplot position
    row = i // 2
    col = i % 2

    # Select the current subplot
    ax = axs[row, col]

    # Filter the data for the current filter type
    filter_data = df_pos_rmse[df_pos_rmse["filter_type"] == filter_type]

    # Get the x-axis positions for the bars
    x_pos = np.arange(len(time_windows))

    # Get the x, y, z RMSE values for the current filter type
    x_rmse = filter_data["x_rmse"]
    y_rmse = filter_data["y_rmse"]
    z_rmse = filter_data["z_rmse"]

    # Plot the x RMSE values as bars
    ax.bar(x_pos, x_rmse, width=0.2, label="x")
    # Plot the y RMSE values as bars with a slight offset
    ax.bar(x_pos + 0.3, y_rmse, width=0.2, label="y")
    # Plot the z RMSE values as bars with a larger offset
    ax.bar(x_pos + 0.6, z_rmse, width=0.2, label="z")

    # Set the x-axis ticks and labels
    ax.set_xticks(x_pos + 0.3)
    # Split the time window into start and end values and format them
    time_window_labels = [
        f"{float(tw.split('-')[0]):.2f}-{float(tw.split('-')[1]):.2f}"
        for tw in time_windows
    ]
    ax.set_xticklabels(time_window_labels)
    ax.set_xlabel("time window (minutes)")
    ax.set_ylabel("rmse (metres)")
    ax.set_title(filter_type)

    # Set the y-axis scale to logarithmic
    ax.set_yscale("log")

    # Set the y-axis limits
    ax.set_ylim([1, max_pos_rmse])

    # Add a legend to the subplot
    ax.legend()

# Adjust the spacing between subplots
plt.tight_layout()

# Show the plot
# plt.show()
fig.savefig("plots/all_pos_rmse_id_" + str(norad_id) + ".pdf")

# Set the figure size and create subplots
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(11, 6))

# Set the font size
plt.rcParams.update({"font.size": 12})

# Calculate the maximum RMSE value across all filter types and time windows
max_vel_rmse = df_vel_rmse[["vx_rmse", "vy_rmse", "vz_rmse"]].max().max() * 1.05

# Loop over the filter types
for i, filter_type in enumerate(filter_types):
    # Calculate the subplot position
    row = i // 2
    col = i % 2

    # Select the current subplot
    ax = axs[row, col]

    # Filter the data for the current filter type
    filter_data = df_vel_rmse[df_vel_rmse["filter_type"] == filter_type]

    # Get the x-axis positions for the bars
    x_vel = np.arange(len(time_windows))

    # Get the x, y, z RMSE values for the current filter type
    vx_rmse = filter_data["vx_rmse"]
    vy_rmse = filter_data["vy_rmse"]
    vz_rmse = filter_data["vz_rmse"]

    # Plot the x RMSE values as bars
    ax.bar(x_vel, vx_rmse, width=0.2, label="vx")
    # Plot the y RMSE values as bars with a slight offset
    ax.bar(x_vel + 0.3, vy_rmse, width=0.2, label="vy")
    # Plot the z RMSE values as bars with a larger offset
    ax.bar(x_vel + 0.6, vz_rmse, width=0.2, label="vz")

    # Set the x-axis ticks and labels
    ax.set_xticks(x_vel + 0.3)
    # Split the time window into start and end values and format them
    time_window_labels = [
        f"{float(tw.split('-')[0]):.2f}-{float(tw.split('-')[1]):.2f}"
        for tw in time_windows
    ]
    ax.set_xticklabels(time_window_labels)
    ax.set_xlabel("time window (minutes)")
    ax.set_ylabel("rmse (m/s)")
    ax.set_title(filter_type)

    # Set the y-axis scale to logarithmic
    ax.set_yscale("log")

    # Set the y-axis limits
    ax.set_ylim([-1, max_vel_rmse])

    # Add a legend to the subplot
    ax.legend()

# Adjust the spacing between subplots
plt.tight_layout()

# Show the plot
# plt.show()
fig.savefig("plots/all_vel_rmse_id_" + str(norad_id) + ".pdf")
