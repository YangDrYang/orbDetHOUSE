import os
import pandas as pd
import matplotlib.pyplot as plt


folder_path = "out/"  # replace with the path to your folder
# Get a sorted list of file names in the folder that start with "ukf_"
trial_file_names = sorted(
    [
        filename
        for filename in os.listdir(folder_path)
        if filename.startswith("meas_") and filename.endswith(".csv")
    ],
)

# Loop through each file
for file_name in trial_file_names:
    file_path = os.path.join(folder_path, file_name)

    # Read the trial CSV file into a pandas dataframe
    trial_df = pd.read_csv(file_path)
    # Add the errors to the empty dataframe
    df["pos_err_x"] = df.get("pos_err_x", 0) + trial_df["EST X1"] - truth_df["x"]


# Read CSV file into a pandas dataframe
df = pd.read_csv("out/" + meas + "_trajectory_error.csv")

# adjust the size of the plot
plt.figure(figsize=(8, 6))
plt.plot(df["time_lapse"], df["pos_err_rms"])

# Add labels and title
plt.xlabel("time_lapse")
plt.ylabel("position errors")
plt.title("position 3D rmse")

# # Show the plot
# plt.show()
plt.savefig("plots/" + filter_type + "_pos_rmse.pdf")

# adjust the size of the plot
plt.figure(figsize=(8, 6))
plt.plot(df["time_lapse"], df["vel_err_rms"])

# Add labels and title
plt.xlabel("time_lapse")
plt.ylabel("velocity errors")
plt.title("velocity 3D rmse")

# # Show the plot
# plt.show()

plt.savefig("plots/" + filter_type + "_vel_rmse.pdf")
