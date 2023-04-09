import pandas as pd
import matplotlib.pyplot as plt
import processing
import seaborn as sns


# filters = ["ukf"]
# filters = ["cut6"]
filters = ["house", "ukf", "cut4", "cut6"]
# filters = ["house", "ukf"]

folder_path = "out_/"

# Create a new figure and axes for logarithmic plot
fig_pos, ax_pos = plt.subplots()

for filter_type in filters:
    processing.process_each_filter(filter_type, folder_path)

    # Read CSV file into a pandas dataframe
    df = pd.read_csv("plots/" + filter_type + "_trajectory_error.csv")

    # Plot the data on the axes
    ax_pos.semilogy(df["time_lapse"], df["pos_err_rms"], label=filter_type)

# Add labels and title
ax_pos.set_xlabel("time_lapse (s)")
ax_pos.set_ylabel("position errors (m)")
ax_pos.set_ylim([-1000, 5000])
ax_pos.set_title("logarithmic 3D position rmse")
ax_pos.legend()
fig_pos.savefig("plots/all_pos_rmse.pdf")

# Create a new figure and axes for violin plot
fig_pos, ax_pos = plt.subplots()
pos_rmse = []
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv("plots/" + filter_type + "_trajectory_error.csv")

    pos_rmse.append(df["pos_err_rms"])

# Create the violin plot
ax_pos.violinplot(pos_rmse, showmeans=False, showmedians=True)
# Add labels and title
ax_pos.set_xlabel("time_lapse (s)")
ax_pos.set_ylabel("position errors (m)")
ax_pos.set_ylim([-1000, 5000])
ax_pos.legend(title=filters, loc="upper right")
ax_pos.set_title("distribution 3D position error")
fig_pos.savefig("plots/all_pos_rmse_violin.pdf")

# Create a new figure and axes for logarithmic plot
fig_vel, ax_vel = plt.subplots()

for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv("plots/" + filter_type + "_trajectory_error.csv")

    # Plot the data on the axes
    ax_vel.plot(df["time_lapse"], df["vel_err_rms"], label=filter_type)

# Add labels and title
ax_vel.set_xlabel("time_lapse (s)")
ax_vel.set_ylabel("velocity errors (m/s)")
ax_vel.set_ylim([-0.5, 5])
ax_vel.set_title("logarithmic 3D velocity rmse")
ax_vel.legend()
fig_vel.savefig("plots/all_vel_rmse.pdf")

# Create a new figure and axes for violin plot
fig_vel, ax_vel = plt.subplots()
vel_rmse = []
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv("plots/" + filter_type + "_trajectory_error.csv")

    vel_rmse.append(df["vel_err_rms"])

# Plot the data on the axes
ax_vel.violinplot(vel_rmse, showmeans=False, showmedians=True)
# Add labels and title
ax_vel.set_xlabel("time_lapse (s)")
ax_vel.set_ylabel("velocity errors (m/s)")
ax_vel.set_ylim([-0.5, 5])
ax_pos.legend(title=filters, loc="upper left")
ax_vel.set_title("distribution of 3D velocity")
fig_vel.savefig("plots/all_vel_rmse_voilin.pdf")
