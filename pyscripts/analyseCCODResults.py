import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import datetime as dt
import scipy.stats as stats
import os

filters = ["srhouse", "house", "srukf", "ukf", "cut4", "cut6"]
# filters = ["ukf"]
# filters = ["srhouse", "house"]
# Directory path
data_folder_path = "out/out_ccdata/"
plot_folder_path = "plots/"
state_type = "mee"
# state_type = "eci"

norad_id = 46984
meas_file = "ccdata/meas_data_id_" + str(norad_id) + ".csv"
stn_file = "ccdata/stn_eci_coordinates.csv"
od_ref_data_file = "refdata/od_ref_id_" + str(norad_id) + ".csv"


# *************** seperate sub windows
meas_df = pd.read_csv(meas_file)
stamp = (meas_df["MJD"] - meas_df["MJD"][0]) * 1440

stamp_diff = np.diff(stamp)
threshold = 10.0
window_starts = np.where(stamp_diff > threshold)[0] + 1

# Add the start and end indices for the first and last time windows
window_starts = np.insert(window_starts, 0, 0)
window_starts = np.append(window_starts, len(stamp))

time_windows = []

for i in range(len(window_starts) - 1):
    start_idx = window_starts[i]
    end_idx = window_starts[i + 1] - 1

    start_time = stamp[start_idx]
    end_time = stamp[end_idx]

    time_windows.append((start_time, end_time))

print("Time Windows:")
for i, (start_time, end_time) in enumerate(time_windows):
    print(f"Time Window {i + 1}: Start = {start_time:.2f}, End = {end_time:.2f}")
    
# # *************** post-residuals

# Create a new figure and axes for plot
num_cols = 2  # Number of filters per row
num_rows = (len(filters) + num_cols - 1) // num_cols  # Calculate the number of rows
fig, axes = plt.subplots(
    ncols=num_cols, nrows=num_rows, figsize=(12, 2 * num_rows), sharex=True
)
# Set the font size
plt.rcParams.update({"font.size": 12})

for i, filter_type in enumerate(filters):
    row = i // num_cols
    col = i % num_cols

    # post measurement residuals
    od_file = (
        data_folder_path
        + filter_type
        + "_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".csv"
    )
    post_res_file = (
        data_folder_path
        + filter_type
        + "_post_res_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".csv"
    )
    post_res_plot_file = (
        plot_folder_path
        + filter_type
        + "_post_res_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".pdf"
    )

    processing.process_post_res_each_filter_ccdata(
        od_file, stn_file, meas_file, post_res_file, post_res_plot_file
    )

    post_res_df = pd.read_csv(post_res_file)
    dates = pd.to_datetime(post_res_df["MJD"] + 2400000.5, unit="D", origin="julian")

    # Plot RA vs MJD
    ax1 = axes[row, col]
    ax2 = ax1.twinx()

    ax1.scatter(
        dates,
        np.degrees(post_res_df["RA"]) * 3600,
        color="blue",
        s=5,
        label="RA Residuals",
    )
    ax1.set_ylabel("RA Residuals \n(arcsec)", color="blue")
    ax1.tick_params(axis="y", labelcolor="blue")
    # ax1.set_ylim(-5e-6, 7.5e-6)
    ax1.set_ylim(-30, 40)
    # ax1.yaxis.set_major_formatter(
    #     mtick.FormatStrFormatter("%.2e")
    # )  # Scientific notation
    if filter_type == "srhouse":
        title = "$w$-HOUSE"
    elif filter_type == "house":
        title = "$\\delta$-HOUSE"
    else:
        title = filter_type.upper()
    # Add filter name as title
    ax1.set_title(title)

    # add 20 mins shift to distinguish ra and dec
    dates = pd.to_datetime(
        post_res_df["MJD"] + 2400000.5 + 20 / 1440, unit="D", origin="julian"
    )
    ax2.scatter(
        dates,
        np.degrees(post_res_df["Dec"]) * 3600,
        color="green",
        s=5,
        label="Dec Residuals",
    )
    ax2.set_ylabel("Dec Residuals \n(arcsec)", color="green")
    ax2.tick_params(axis="y", labelcolor="green")
    # ax2.set_ylim(-1e-5, 1.0e-5)
    ax2.set_ylim(-25, 40)
    # ax2.yaxis.set_major_formatter(
    #     mtick.FormatStrFormatter("%.2e")
    # )  # Scientific notation

    # Format x-axis as hours
    hours = mdates.HourLocator(interval=6)
    hour_format = mdates.DateFormatter("%H:%M")
    ax2.xaxis.set_major_locator(hours)
    ax2.xaxis.set_major_formatter(hour_format)

    # Show only the date for the first epoch of each day
    dates_first_epoch = np.unique([date.date() for date in dates])
    for date in dates_first_epoch:
        ax1.axvline(date, color="gray", linestyle="--", alpha=0.5)
        ax2.axvline(date, color="gray", linestyle="--", alpha=0.5)

    # Show only 1 hour before and 1 hour after the data on x-axis
    start_time = min(dates) - dt.timedelta(hours=1)
    end_time = max(dates) + dt.timedelta(hours=1)
    ax2.set_xlim(start_time, end_time)

    # # Add legend
    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2, loc="upper left")
    
    # Convert RA and Dec from degrees to arcseconds
    ra_arcsec = np.degrees(post_res_df["RA"]) * 3600
    dec_arcsec = np.degrees(post_res_df["Dec"]) * 3600

    # Add these to your DataFrame
    post_res_df["RA_arcsec"] = ra_arcsec
    post_res_df["Dec_arcsec"] = dec_arcsec

    # Convert MJD to minutes from start
    post_res_df["MJD_minutes"] = (post_res_df["MJD"] - post_res_df["MJD"][0]) * 1440

    # Loop over time windows
    for i, (start_time, end_time) in enumerate(time_windows):
        # Filter DataFrame for current time window
        window_df = post_res_df[(post_res_df["MJD_minutes"] >= start_time) & (post_res_df["MJD_minutes"] <= end_time)]

        # Calculate the square of residuals
        ra_squared = np.square(window_df["RA_arcsec"])
        dec_squared = np.square(window_df["Dec_arcsec"])

        # Calculate the mean of squared residuals
        ra_mean_squared = np.mean(ra_squared)
        dec_mean_squared = np.mean(dec_squared)

        # Calculate the RMS of residuals
        ra_rms = np.sqrt(ra_mean_squared)
        dec_rms = np.sqrt(dec_mean_squared)

        # Print the RMS of residuals
        print(f"Time Window {i + 1}: Start = {start_time:.2f}, End = {end_time:.2f}")
        print(f"RMS of RA residuals for {filter_type}: {ra_rms}")
        print(f"RMS of Dec residuals for {filter_type}: {dec_rms}")


# Rotate x-axis tick labels and adjust label spacing
fig.autofmt_xdate(rotation=45)
plt.tight_layout()

post_res_plot_all_file = (
    plot_folder_path + "all_post_res_id_" + str(norad_id) + "_" + state_type + ".pdf"
)
# Save the figure
plt.savefig(post_res_plot_all_file)


filters = ["house", "srukf", "ukf", "cut4", "cut6", "srhouse"]
# *************** absolute error plot for position components

# Create an empty DataFrame to store RMS values
df_pos_rmse = pd.DataFrame(
    columns=["filter_type", "time_window", "x_rmse", "y_rmse", "z_rmse", "pos_3d_rmse"]
)

# Plot the data with the identified time windows
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))

# Set the font size
plt.rcParams.update({"font.size": 12})

# Create a color mapping dictionary
color_map = {
    "srhouse": "green",
    "house": "black",
    "srukf": "pink",
    "ukf": "red",
    "cut4": "blue",
    "cut6": "orange",
}

for i, filter_type in enumerate(filters):
    processing.process_rmse_each_filter_ccdata(
        filter_type, data_folder_path, norad_id, od_ref_data_file, state_type
    )
    # Read the data into a pandas dataframe
    err_df = pd.read_csv(
        data_folder_path
        + filter_type
        + "_err_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".csv"
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
        rmse = np.sqrt(
            np.mean(x_errors**2) + np.mean(y_errors**2) + np.mean(z_errors**2)
        )

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
                        "pos_3d_rmse": [rmse],
                    }
                ),
            ],
            ignore_index=True,
        )

        if filter_type == "srhouse":
            label = "$w$-HOUSE"
        elif filter_type == "house":
            label = "$\\delta$-HOUSE"
        else:
            label = filter_type.upper()
        sub_ax = ax[0, j]
        if j == 0:
            sub_ax.scatter(
                stamp[start_idx],
                np.abs(err_df["pos_err_x"][start_idx]),
                color=color,
                label=label,
            )
        else:
            sub_ax.scatter(
                stamp[start_idx],
                np.abs(err_df["pos_err_x"][start_idx]),
                color=color,
            )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("x absolute error (m)")
        sub_ax.set_yscale("log")
        if j == 0:
            sub_ax.legend()

    # Plot absolute y errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[1, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["pos_err_y"][start_idx]),
            color=color,
            # label=filter_type.upper(),
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("y absolute error (m)")
        sub_ax.set_yscale("log")
        # sub_ax.legend()

    # Plot absolute z errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[2, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["pos_err_z"][start_idx]),
            color=color,
            # label=filter_type.upper(),
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("z absolute error (m)")
        sub_ax.set_yscale("log")
        # sub_ax.legend()

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
fig.savefig("plots/all_pos_abserr_id_" + str(norad_id) + "_" + state_type + ".pdf")
# Print the DataFrame with RMS values
# Set the number of decimal places to 2 for printing
pd.options.display.float_format = "{:.3f}".format
print(df_pos_rmse)


# *************** absolute error plot for velocity components

# Create an empty DataFrame to store RMS values
df_vel_rmse = pd.DataFrame(
    columns=[
        "filter_type",
        "time_window",
        "vx_rmse",
        "vy_rmse",
        "vz_rmse",
        "vel_3d_rmse",
    ]
)

# Plot the data with the identified time windows
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))

# Set the font size
plt.rcParams.update({"font.size": 12})

for i, filter_type in enumerate(filters):
    # Read the data into a pandas dataframe
    err_df = pd.read_csv(
        data_folder_path
        + filter_type
        + "_err_id_"
        + str(norad_id)
        + "_"
        + state_type
        + ".csv"
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
        rmse = np.sqrt(
            np.mean(vx_errors**2) + np.mean(vy_errors**2) + np.mean(vz_errors**2)
        )

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
                        "vel_3d_rmse": [rmse],
                    }
                ),
            ],
            ignore_index=True,
        )

        if filter_type == "srhouse":
            label = "$w$-HOUSE"
        elif filter_type == "house":
            label = "$\\delta$-HOUSE"
        else:
            label = filter_type.upper()
        sub_ax = ax[0, j]
        if j == 0:
            sub_ax.scatter(
                stamp[start_idx],
                np.abs(err_df["vel_err_x"][start_idx]),
                color=color,
                label=label,
            )
        else:
            sub_ax.scatter(
                stamp[start_idx],
                np.abs(err_df["vel_err_x"][start_idx]),
                color=color,
                # label=filter_type.upper(),
            )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vx absolute error (m)")
        sub_ax.set_yscale("log")
        if j == 0:
            sub_ax.legend()

    # Plot absolute y errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[1, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["vel_err_y"][start_idx]),
            color=color,
            # label=filter_type.upper(),
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vy absolute error (m)")
        sub_ax.set_yscale("log")
        # sub_ax.legend()

    # Plot absolute z errors in desired sub-time windows
    for j, (start_time, end_time) in enumerate(time_windows):
        sub_ax = ax[2, j]
        start_idx = (stamp >= start_time) & (stamp <= end_time)
        sub_ax.scatter(
            stamp[start_idx],
            np.abs(err_df["vel_err_z"][start_idx]),
            color=color,
            # label=filter_type.upper(),
        )
        sub_ax.set_xlabel("time (minutes)")
        sub_ax.set_ylabel("vz absolute error (m)")
        sub_ax.set_yscale("log")
        # sub_ax.legend()

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
fig.savefig("plots/all_vel_abserr_id_" + str(norad_id) + "_" + state_type + ".pdf")
# Print the DataFrame with RMS values
print(df_vel_rmse)

# # *************** rmse bar chart for position components

# # Define the filter types and time windows
# time_windows = [
#     "0.0-2.1686866565141827",
#     "1462.6672116550617-1464.8276216641534",
#     "2926.092338307062-2927.989939986728",
# ]

# # Set the figure size and create subplots
# fig, axs = plt.subplots(nrows=len(filters) // 2, ncols=2, figsize=(11, 6))

# # Set the font size
# plt.rcParams.update({"font.size": 12})

# # Calculate the maximum RMSE value across all filter types and time windows
# max_pos_rmse = df_pos_rmse[["x_rmse", "y_rmse", "z_rmse"]].max().max() * 1.05

# # Loop over the filter types
# for i, filter_type in enumerate(filters):
#     # Calculate the subplot position
#     row = i // 2
#     col = i % 2

#     # Select the current subplot
#     ax = axs[row, col]

#     # Filter the data for the current filter type
#     filter_data = df_pos_rmse[df_pos_rmse["filter_type"] == filter_type]

#     # Get the x-axis positions for the bars
#     x_pos = np.arange(len(time_windows))
#     # print("x_pos:   ", x_pos)

#     # Get the x, y, z RMSE values for the current filter type
#     x_rmse = filter_data["x_rmse"]
#     y_rmse = filter_data["y_rmse"]
#     z_rmse = filter_data["z_rmse"]
#     # print("x_rmse:   ", x_rmse)

#     # Plot the x RMSE values as bars
#     ax.bar(x_pos, x_rmse, width=0.2, label="x")
#     # Plot the y RMSE values as bars with a slight offset
#     ax.bar(x_pos + 0.3, y_rmse, width=0.2, label="y")
#     # Plot the z RMSE values as bars with a larger offset
#     ax.bar(x_pos + 0.6, z_rmse, width=0.2, label="z")

#     # Set the x-axis ticks and labels
#     ax.set_xticks(x_pos + 0.3)
#     # Split the time window into start and end values and format them
#     time_window_labels = [
#         f"{float(tw.split('-')[0]):.2f}-{float(tw.split('-')[1]):.2f}"
#         for tw in time_windows
#     ]
#     ax.set_xticklabels(time_window_labels)
#     ax.set_xlabel("time window (minutes)")
#     ax.set_ylabel("rmse (metres)")
#     if i == 0:
#         ax.set_title(filter_type.upper())

#     # Set the y-axis scale to logarithmic
#     ax.set_yscale("log")

#     # Set the y-axis limits
#     ax.set_ylim([1, max_pos_rmse])

#     # Add a legend to the subplot
#     ax.legend(loc="upper left")

# # Adjust the spacing between subplots
# plt.tight_layout()

# # Show the plot
# # plt.show()
# fig.savefig("plots/all_pos_rmse_id_" + str(norad_id) + "_" + state_type + ".pdf")


# # *************** rmse bar chart for velocity components
# # Set the figure size and create subplots
# fig, axs = plt.subplots(nrows=len(filters) // 2, ncols=2, figsize=(12, 6))

# # Set the font size
# plt.rcParams.update({"font.size": 12})

# # Calculate the maximum RMSE value across all filter types and time windows
# max_vel_rmse = df_vel_rmse[["vx_rmse", "vy_rmse", "vz_rmse"]].max().max() * 1.05

# # Loop over the filter types
# for i, filter_type in enumerate(filters):
#     # Calculate the subplot position
#     row = i // 2
#     col = i % 2

#     # Select the current subplot
#     ax = axs[row, col]

#     # Filter the data for the current filter type
#     filter_data = df_vel_rmse[df_vel_rmse["filter_type"] == filter_type]

#     # Get the x-axis positions for the bars
#     x_vel = np.arange(len(time_windows))

#     # Get the x, y, z RMSE values for the current filter type
#     vx_rmse = filter_data["vx_rmse"]
#     vy_rmse = filter_data["vy_rmse"]
#     vz_rmse = filter_data["vz_rmse"]

#     # Plot the x RMSE values as bars
#     ax.bar(x_vel, vx_rmse, width=0.2, label="vx")
#     # Plot the y RMSE values as bars with a slight offset
#     ax.bar(x_vel + 0.3, vy_rmse, width=0.2, label="vy")
#     # Plot the z RMSE values as bars with a larger offset
#     ax.bar(x_vel + 0.6, vz_rmse, width=0.2, label="vz")

#     # Set the x-axis ticks and labels
#     ax.set_xticks(x_vel + 0.3)
#     # Split the time window into start and end values and format them
#     time_window_labels = [
#         f"{float(tw.split('-')[0]):.2f}-{float(tw.split('-')[1]):.2f}"
#         for tw in time_windows
#     ]
#     ax.set_xticklabels(time_window_labels)
#     ax.set_xlabel("time window (minutes)")
#     ax.set_ylabel("rmse (m/s)")
#     if i == 0:
#         ax.set_title(filter_type.upper())

#     # Set the y-axis scale to logarithmic
#     ax.set_yscale("log")

#     # Set the y-axis limits
#     ax.set_ylim([0.01, max_vel_rmse])

#     # Add a legend to the subplot
#     ax.legend(loc="upper left")

# # Adjust the spacing between subplots
# plt.tight_layout()

# # Show the plot
# # plt.show()
# fig.savefig("plots/all_vel_rmse_id_" + str(norad_id) + "_" + state_type + ".pdf")

# # *************** normalised error square
# # flag to determine to process the NES or use existing files for plots
# filters = ["srhouse", "house", "srukf", "ukf", "cut4", "cut6"]
# flag_proc = 1

# # Set the alpha value
# alpha = 0.05
# # Calculate the degrees of freedom
# dof = 6 - 1
# # Calculate the chi-square value at the alpha/2 probability level
# chi2_lower = stats.chi2.ppf(alpha / 2, dof)
# # Calculate the chi-square value at the 1-alpha/2 probability level
# chi2_upper = stats.chi2.ppf(1 - alpha / 2, dof)

# # Create a new figure and axes for logarithmic plot
# fig, ax = plt.subplots(ncols=1, figsize=(5, 4))

# for filter_type in filters:
#     if flag_proc:
#         processing.process_nes_each_filter_ccdata(
#             filter_type, data_folder_path, norad_id, od_ref_data_file, state_type
#         )

#     # Read CSV file into a pandas dataframe
#     df = pd.read_csv(
#         data_folder_path
#         + "nes_"
#         + filter_type
#         + "_id_"
#         + str(norad_id)
#         + "_"
#         + state_type
#         + ".csv"
#     )

#     # Plot the data on the axes
#     ax.scatter(df["tSec"], df["NES"], s=3, label=filter_type)

#     if filter_type == "cut6":
#         # Add a horizontal line for the upper bound
#         ax.axhline(y=chi2_upper, color="r", linestyle="--", label="upper bound")
#         # Add a horizontal line for the lower bound
#         ax.axhline(y=chi2_lower, color="k", linestyle="--", label="lower bound")
#     else:
#         # Add a horizontal line for the upper bound
#         ax.axhline(y=chi2_upper, color="r", linestyle="--")
#         # Add a horizontal line for the lower bound
#         ax.axhline(y=chi2_lower, color="k", linestyle="--")

# # Add labels and title
# ax.set_xlabel("time elapsed (s)")
# ax.set_ylabel("normalised error square (m)")
# # ax.set_ylim([-1, 15])
# # ax.set_title("nes")
# ax.legend()

# plt.tight_layout()
# fig.savefig(
#     plot_folder_path + "all_nes_id_" + str(norad_id) + "_" + state_type + ".pdf"
# )

# *************** plot run time for all filters

# List of filters
filters = ["srhouse", "house", "srukf", "ukf", "cut4", "cut6"]

# Create an empty DataFrame to store the run times
df_run_times = pd.DataFrame()

# Loop over each filter
for filter_type in filters:
    # Construct the file name for the current filter
    file_name = file_name = f"run_times_{filter_type}_id_{norad_id}_{state_type}.csv"

    # Read the CSV file into a DataFrame
    file_path = os.path.join(data_folder_path, file_name)
    filter_df = pd.read_csv(file_path)

    if filter_type == "srhouse":
        label = "$w$-HOUSE"
    elif filter_type == "house":
        label = "$\\delta$-HOUSE"
    else:
        label = filter_type.upper()
    # Extract the run times column and add it to the main DataFrame
    df_run_times[label] = filter_df[filter_type]

print(df_run_times)
# Create a plot
plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(figsize=(8, 4))
df_run_times.plot(kind="barh", ax=ax)

ax.set_ylabel("Filter")  # Set ylabel to represent filters
ax.set_yticklabels([])  # Remove y-axis tick labels
ax.set_xlabel("Run Time (second)")  # Set xlabel to represent run time
ax.set_title("Run Time for Different Filters")

plt.tight_layout()
plt.savefig(
    plot_folder_path + "run_times_id_" + str(norad_id) + "_" + state_type + ".pdf"
)
