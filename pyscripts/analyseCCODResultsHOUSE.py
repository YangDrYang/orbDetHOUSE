import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
import scipy.stats as stats

filters = ["house"]
# Directory path
out_folder_path = "out/out_ccdata/"
plot_folder_path = "plots/"

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


# *************** absolute error plot for position components

trial_total_num = 100
# Create an empty DataFrame to store RMS values
df_pos_rmse = pd.DataFrame(
    columns=["trial_num", "time_window", "x_rmse", "y_rmse", "z_rmse", "3d_rmse"]
)

# Plot the data with the identified time windows
# fig, ax = plt.subplots(nrows=3, ncols=len(time_windows), figsize=(15, 10))
fig, ax = plt.subplots(nrows=1, ncols=len(time_windows), figsize=(15, 5))

for j, (start_time, end_time) in enumerate(time_windows):
    start_idx = (stamp >= start_time) & (stamp <= end_time)
    sub_ax = ax[j]  # Adjust subplot indexing
    for i in range(trial_total_num):
        trial_no = i + 1
        processing.process_rmse_HOUSE_ccdata(
            trial_no, out_folder_path, norad_id, od_ref_data_file
        )
        # Read the data into a pandas dataframe
        err_df = pd.read_csv(
            out_folder_path
            + "house_err_id_"
            + str(norad_id)
            + "_"
            + str(trial_no)
            + ".csv"
        )

        # Extract the absolute errors within the current time window
        x_errors = np.abs(err_df["pos_err_x"][start_idx])
        y_errors = np.abs(err_df["pos_err_y"][start_idx])
        z_errors = np.abs(err_df["pos_err_z"][start_idx])

        # Calculate RMS for each error component within the sub-time window
        x_rmse = np.sqrt(np.mean(x_errors**2))
        y_rmse = np.sqrt(np.mean(y_errors**2))
        z_rmse = np.sqrt(np.mean(z_errors**2))
        rmse = np.sqrt(np.mean(x_errors**2 + y_errors**2 + z_errors**2))

        # Store the RMS values in the DataFrame
        df_pos_rmse = pd.concat(
            [
                df_pos_rmse,
                pd.DataFrame(
                    {
                        "trial_num": [trial_no],
                        "time_window": [f"{start_time}-{end_time}"],
                        "x_rmse": [x_rmse],
                        "y_rmse": [y_rmse],
                        "z_rmse": [z_rmse],
                        "3d_rmse": [rmse],
                    }
                ),
            ],
            ignore_index=True,
        )

        # sub_ax_x = ax[0, j]
        # sub_ax_y = ax[1, j]
        # sub_ax_z = ax[2, j]

        # # Plot x absolute errors
        # sub_ax_x.scatter(
        #     stamp[start_idx],
        #     x_errors,
        #     # label=f"Trial {i}",
        #     marker=".",
        # )
        # sub_ax_x.set_xlabel("time (minutes)")
        # sub_ax_x.set_ylabel("x absolute error (m)")
        # sub_ax_x.set_yscale("log")
        # # sub_ax_x.legend()

        # # Plot y absolute errors
        # sub_ax_y.scatter(
        #     stamp[start_idx],
        #     y_errors,
        #     # label=f"Trial {i}",
        #     marker=".",
        # )
        # sub_ax_y.set_xlabel("time (minutes)")
        # sub_ax_y.set_ylabel("y absolute error (m)")
        # sub_ax_y.set_yscale("log")
        # # sub_ax_y.legend()

        # # Plot z absolute errors
        # sub_ax_z.scatter(
        #     stamp[start_idx],
        #     z_errors,
        #     # label=f"Trial {i}",
        #     marker=".",
        # )
        # sub_ax_z.set_xlabel("time (minutes)")
        # sub_ax_z.set_ylabel("z absolute error (m)")
        # sub_ax_z.set_yscale("log")
        # # sub_ax_z.legend()

    # Plot the 3D RMSE for all trials in a line plot
    sub_ax.plot(
        range(trial_total_num),
        df_pos_rmse[df_pos_rmse["time_window"] == f"{start_time}-{end_time}"][
            "3d_rmse"
        ],
        linestyle="-",
        linewidth=2,
        marker="o",
        markersize=4,
    )

    sub_ax.set_xlabel("trial number")
    sub_ax.set_ylabel("3d rmse")
    sub_ax.set_title(f"time window: {start_time:.3f} min - {end_time:.3f} min ")

# # Set x-axis limits and x-ticks for each subplot to show all three time windows
# for j, (start_time, end_time) in enumerate(time_windows):
#     sub_ax_x = ax[0, j]
#     sub_ax_y = ax[1, j]
#     sub_ax_z = ax[2, j]
#     sub_ax_x.set_xlim(start_time, end_time)
#     sub_ax_x.set_xticks([start_time, end_time])
#     sub_ax_y.set_xlim(start_time, end_time)
#     sub_ax_y.set_xticks([start_time, end_time])
#     sub_ax_z.set_xlim(start_time, end_time)
#     sub_ax_z.set_xticks([start_time, end_time])

# Set font size for tick labels
for sub_ax in ax.flat:
    sub_ax.tick_params(axis="both", which="both", labelsize=12)

plt.tight_layout()
fig.savefig("plots/all_pos_abserr_id_" + str(norad_id) + "_all_thetas.pdf")

# Print the DataFrame with RMS values
print(df_pos_rmse)
df_pos_rmse.to_csv(
    out_folder_path + "house_err_id_" + str(norad_id) + "_all_thetas.csv",
    index=False,
)


# # *************** absolute error plot for velocity components

# # Create an empty DataFrame to store RMS values
# df_vel_rmse = pd.DataFrame(
#     columns=["filter_type", "time_window", "vx_rmse", "vy_rmse", "vz_rmse"]
# )

# # Plot the data with the identified time windows
# fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))

# for i, filter_type in enumerate(filters):
#     # Read the data into a pandas dataframe
#     err_df = pd.read_csv(
#         out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
#     )

#     # Get the color for the current filter type from the color mapping dictionary
#     color = color_map.get(
#         filter_type, "black"
#     )  # Use black as the default color if not found in the mapping

#     # Plot absolute x errors in desired sub-time windows
#     for j, (start_time, end_time) in enumerate(time_windows):
#         start_idx = (stamp >= start_time) & (stamp <= end_time)
#         vx_errors = err_df["vel_err_x"][start_idx]
#         vy_errors = err_df["vel_err_y"][start_idx]
#         vz_errors = err_df["vel_err_z"][start_idx]

#         # Calculate RMS for each error component within the sub-time window
#         vx_rmse = np.sqrt(np.mean(vx_errors**2))
#         vy_rmse = np.sqrt(np.mean(vy_errors**2))
#         vz_rmse = np.sqrt(np.mean(vz_errors**2))

#         # Store the RMS values in the DataFrame
#         df_vel_rmse = pd.concat(
#             [
#                 df_vel_rmse,
#                 pd.DataFrame(
#                     {
#                         "filter_type": [filter_type],
#                         "time_window": [f"{start_time}-{end_time}"],
#                         "vx_rmse": [vx_rmse],
#                         "vy_rmse": [vy_rmse],
#                         "vz_rmse": [vz_rmse],
#                     }
#                 ),
#             ],
#             ignore_index=True,
#         )

#         sub_ax = ax[0, j]
#         sub_ax.scatter(
#             stamp[start_idx],
#             np.abs(err_df["vel_err_x"][start_idx]),
#             color=color,
#             label=filter_type,
#         )
#         sub_ax.set_xlabel("time (minutes)")
#         sub_ax.set_ylabel("vx absolute error (m)")
#         sub_ax.set_yscale("log")
#         sub_ax.legend()

#     # Plot absolute y errors in desired sub-time windows
#     for j, (start_time, end_time) in enumerate(time_windows):
#         sub_ax = ax[1, j]
#         start_idx = (stamp >= start_time) & (stamp <= end_time)
#         sub_ax.scatter(
#             stamp[start_idx],
#             np.abs(err_df["vel_err_y"][start_idx]),
#             color=color,
#             label=filter_type,
#         )
#         sub_ax.set_xlabel("time (minutes)")
#         sub_ax.set_ylabel("vy absolute error (m)")
#         sub_ax.set_yscale("log")
#         sub_ax.legend()

#     # Plot absolute z errors in desired sub-time windows
#     for j, (start_time, end_time) in enumerate(time_windows):
#         sub_ax = ax[2, j]
#         start_idx = (stamp >= start_time) & (stamp <= end_time)
#         sub_ax.scatter(
#             stamp[start_idx],
#             np.abs(err_df["vel_err_z"][start_idx]),
#             color=color,
#             label=filter_type,
#         )
#         sub_ax.set_xlabel("time (minutes)")
#         sub_ax.set_ylabel("vz absolute error (m)")
#         sub_ax.set_yscale("log")
#         sub_ax.legend()

# # Set x-axis limits and x-ticks for each subplot to show all three time windows
# for i in range(3):
#     for j, (start_time, end_time) in enumerate(time_windows):
#         sub_ax = ax[i, j]
#         sub_ax.set_xlim(start_time, end_time)
#         sub_ax.set_xticks([start_time, end_time])

# # Set font size for tick labels
# for sub_ax in ax.flat:
#     sub_ax.tick_params(axis="both", which="both", labelsize=12)
# plt.tight_layout()
# # plt.show()
# fig.savefig("plots/all_vel_abserr_id_" + str(norad_id) + ".pdf")
# # Print the DataFrame with RMS values
# print(df_vel_rmse)
