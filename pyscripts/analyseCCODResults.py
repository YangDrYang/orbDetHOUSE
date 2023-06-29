import pandas as pd
import matplotlib.pyplot as plt
import processing

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
# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(ncols=2, figsize=(10, 5))


for filter_type in filters:
    # post measurement residuals
    od_file = "out_ccdata/" + filter_type + "_id_" + str(norad_id) + ".csv"
    # processing.process_post_res_each_filter(od_file, stn_file, meas_file, post_res_file)

    processing.process_rmse_each_filter_ccdata(
        filter_type, out_folder_path, norad_id, od_ref_data_file
    )

    # Read CSV file into a pandas dataframe
    err_df = pd.read_csv(
        out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
    )

    stamp = (err_df["mjd"] - err_df["mjd"][0]) * 1440

    # Plot the data on the axes
    ax[0].plot(stamp, err_df["pos_err_rms"], label=filter_type)

    # Plot the data on the axes
    ax[1].plot(stamp, err_df["vel_err_rms"], label=filter_type)

# Add labels and title
ax[0].set_xlabel("time elapsed (s)")
ax[0].set_ylabel("position errors (m)")
# ax[0].set_ylim([-10000, 10000])
ax[0].set_title("3D position rmse")
ax[0].legend()

# Add labels and title
ax[1].set_xlabel("time elapsed (s)")
ax[1].set_ylabel("velocity errors (m/s)")
# ax[1].set_ylim([-1, 1])
ax[1].set_title("3D velocity rmse")
ax[1].legend()

plt.tight_layout()
fig.savefig("plots/all_rmse_id_" + str(norad_id) + ".pdf")

# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(nrows=3, figsize=(10, 5))
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    err_df = pd.read_csv(
        out_folder_path + filter_type + "_err_id_" + str(norad_id) + ".csv"
    )
    stamp = (err_df["mjd"] - err_df["mjd"][0]) * 1440
    # Plot the data on the axes
    # ax[0].plot(err_df["mjd"].head(70))
    # ax[1].plot(err_df["pos_err_x"].head(70))
    # ax[2].plot(err_df["pos_err_y"].head(70))
    # ax[3].plot(err_df["pos_err_z"].head(70))
    # ax[0].scatter(err_df["mjd"])
    ax[0].plot(stamp, err_df["pos_err_x"])
    ax[1].plot(stamp, err_df["pos_err_y"])
    ax[2].plot(stamp, err_df["pos_err_z"])

# Add labels and title
# ax[0].set_ylabel("mjd (day)")
# Add labels and title
ax[0].set_ylabel("x errors (m)")
# ax[0].set_ylim([-10000, 10000])
ax[0].legend()
# Add labels and title
ax[1].set_ylabel("y errors (m)")
# ax[1].set_ylim([-10000, 10000])
ax[1].legend()
# Add labels and title
ax[2].set_ylabel("z errors (m)")
# ax[2].set_ylim([-10000, 10000])
ax[2].legend()

plt.tight_layout()
plt.show()
