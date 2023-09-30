import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
import scipy.stats as stats
import seaborn as sns

filters = ["house", "srhouse"]
filter_type = "house"
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

# Create an empty DataFrame to store RMS values
df_pos_rmse = pd.DataFrame(
    columns=["trial_num", "time_window", "x_rmse", "y_rmse", "z_rmse", "3d_rmse"]
)

# Plot the data with the identified time windows
# fig, ax = plt.subplots(nrows=3, ncols=len(time_windows), figsize=(15, 10))
fig, ax = plt.subplots(nrows=len(time_windows), ncols=1, figsize=(8, 6))

for j, (start_time, end_time) in enumerate(time_windows):
    start_idx = (stamp >= start_time) & (stamp <= end_time)
    sub_ax = ax[j]  # Adjust subplot indexing
    trial_total_num = 100
    for i in range(trial_total_num):
        trial_no = i + 1
        processing.process_rmse_HOUSE_ccdata(
            filter_type, trial_no, out_folder_path, norad_id, od_ref_data_file
        )
        try:
            # Read the data into a pandas dataframe
            err_df = pd.read_csv(
                out_folder_path
                + "house_err_id_"
                + str(norad_id)
                + "_"
                + str(trial_no)
                + ".csv"
            )
        except FileNotFoundError:
            trial_total_num = trial_total_num - 1
            continue

        # Extract the absolute errors within the current time window
        x_errors = np.abs(err_df["pos_err_x"][start_idx])
        y_errors = np.abs(err_df["pos_err_y"][start_idx])
        z_errors = np.abs(err_df["pos_err_z"][start_idx])

        # Calculate RMS for each error component within the sub-time window
        x_rmse = np.sqrt(np.mean(x_errors**2))
        y_rmse = np.sqrt(np.mean(y_errors**2))
        z_rmse = np.sqrt(np.mean(z_errors**2))
        rmse = np.sqrt(np.mean(x_errors**2 + y_errors**2 + z_errors**2))

        # Format start_time and end_time with two decimal places
        formatted_start_time = f"{start_time:.2f}"
        formatted_end_time = f"{end_time:.2f}"
        # Create the "time_window" column with formatted times
        time_window = f"{formatted_start_time}-{formatted_end_time}"

        # Store the RMS values in the DataFrame
        df_pos_rmse = pd.concat(
            [
                df_pos_rmse,
                pd.DataFrame(
                    {
                        "trial_num": [trial_no],
                        # "time_window": [f"{start_time}-{end_time}"],
                        "time_window": [time_window],
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
        df_pos_rmse[df_pos_rmse["time_window"] == f"{start_time:.2f}-{end_time:.2f}"][
            "3d_rmse"
        ],
        linestyle="-",
        linewidth=2,
        marker="o",
        markersize=4,
    )
    # Set the y-axis scale to logarithmic
    sub_ax.set_yscale("log")
    sub_ax.set_xlabel("trial number")
    sub_ax.set_ylabel("3d rmse")
    sub_ax.set_title(f"time window: {start_time:.3f} min - {end_time:.3f} min ")


# Set font size for tick labels
for sub_ax in ax.flat:
    sub_ax.tick_params(axis="both", which="both", labelsize=12)

plt.tight_layout()
file_name = (
    "plots/" + filter_type + "_3d_rmse_vs_theta_id_" + str(norad_id) + "_plot.pdf"
)
fig.savefig(file_name)

# Print the DataFrame with RMS values
print(df_pos_rmse)
file_name = (
    out_folder_path + filter_type + "_err_id_" + str(norad_id) + "_all_thetas.csv"
)
df_pos_rmse.to_csv(
    file_name,
    index=False,
)

# Extract trial_nums
trial_nums = df_pos_rmse["trial_num"].unique().tolist()

# Extract time_windows
time_windows = df_pos_rmse["time_window"].unique().tolist()

# Initialize an empty dictionary to store rmse_values
rmse_values = {}

# Iterate over each row of the DataFrame
for _, row in df_pos_rmse.iterrows():
    time_window = row["time_window"]
    rmse_value = row["3d_rmse"]

    # Check if time_window exists as a key in rmse_values dictionary
    if time_window in rmse_values:
        rmse_values[time_window].append(rmse_value)
    else:
        rmse_values[time_window] = [rmse_value]

# Set up the figure and axis
fig, ax = plt.subplots()

# Set the width of each bar
bar_width = 0.2

# Set the x positions of the bars
x_pos = np.arange(len(trial_nums))

# Plot each group of bars
for i, time_window in enumerate(time_windows):
    if time_window in rmse_values:
        ax.bar(
            x_pos + (i * bar_width),
            rmse_values[time_window],
            bar_width,
            label=f"{time_window} min",
        )

# Set the x-axis ticks and labels
ax.set_xticks(x_pos + (len(time_windows) / 2) * bar_width)
ax.set_xticklabels(trial_nums)
ax.set_xlabel("Trial Number")

# Set the y-axis label
ax.set_ylabel("3D RMSE (m)")

# Set the y-axis scale to logarithmic
ax.set_yscale("log")

# Set the chart title
ax.set_title("3D RMSE of $\delta$-HOUSE (m)")

# Add a legend
ax.legend()

# Save the figure as an image file
file_name = (
    "plots/" + filter_type + "_3d_rmse_vs_theta_id_" + str(norad_id) + "_bar.pdf"
)
plt.savefig(file_name)  # Change the filename and extension as needed
