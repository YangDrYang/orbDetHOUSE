import pandas as pd
import matplotlib.pyplot as plt
import processing

filters = ["house", "ukf", "cut4", "cut6"]
# Directory path
folder_path = "out_ccdata"

norad_id = 46984
meas_file = "ccdata/meas_data_id_" + str(norad_id) + ".csv"
stn_file = "ccdata/stn_data_id_" + str(norad_id) + ".csv"
post_res_file = "out_ccdata/post_res_id_" + str(norad_id) + ".csv"
# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(ncols=2, figsize=(10, 5))

for filter_type in filters:
    od_file = "out_ccdata/" + filter_type + ".csv"
    processing.process_post_res_each_filter(od_file, stn_file, meas_file, post_res_file)

    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "trajectory_error_" + filter_type + ".csv")

    # Plot the data on the axes
    ax[0].semilogy(df["tSec"], df["pos_err_rms"], label=filter_type)

    # Plot the data on the axes
    ax[1].semilogy(df["tSec"], df["vel_err_rms"], label=filter_type)

# Add labels and title
ax[0].set_xlabel("time elapsed (s)")
ax[0].set_ylabel("position errors (m)")
# ax_pos.set_ylim([-1000, 5000])
ax[0].set_title("3D position rmse")
ax[0].legend()

# Add labels and title
ax[1].set_xlabel("time elapsed (s)")
ax[1].set_ylabel("velocity errors (m/s)")
# ax[1].set_ylim([-0.5, 1])
ax[1].set_title("3D velocity rmse")
ax[1].legend()

plt.tight_layout()
fig.savefig("plots/all_rmse.pdf")

# post measurement residuals
