import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


folder_path = "out_/"  # replace with the path to your folder
# Get a sorted list of file names in the folder that start with "ukf_"
# trial_file_names = sorted(
#     [
#         filename
#         for filename in os.listdir(folder_path)
#         if filename.startswith("meas_") and filename.endswith(".csv")
#     ],
# )
df = pd.read_csv(folder_path + "measurement_truth.csv")

NO_MEASUREMENT = 99999999
df = df.replace(NO_MEASUREMENT, np.nan)

# Extract the required columns
tSec = df["tSec"]
ra = df["ra"]
dec = df["dec"]
range = df["range"]
range_rate = df["range_rate"]


# Create the subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Plot ra vs tSec
axs[0, 0].scatter(tSec, ra / np.pi * 180, s=2)
axs[0, 0].set_xlim(0, tSec[len(tSec) - 1])
# axs[0,0].set_xlabel("tSec (seconds)")
axs[0, 0].set_ylabel("ra (degrees)")

# Plot dec vs tSec
axs[0, 1].scatter(tSec, dec / np.pi * 180, s=2)
axs[0, 1].set_xlim(0, tSec[len(tSec) - 1])
# axs[0,1].set_xlabel("tSec")
axs[0, 1].set_ylabel("dec (degrees)")

# Plot range vs tSec
axs[1, 0].scatter(tSec, range / 1000, s=2)
axs[1, 0].set_xlim(0, tSec[len(tSec) - 1])
axs[1, 0].set_xlabel("tSec (seconds)")
axs[1, 0].set_ylabel("range (km)")

# Plot range_rate vs tSec
axs[1, 1].scatter(tSec, range_rate / 1000, s=2)
axs[1, 1].set_xlim(0, tSec[len(tSec) - 1])
axs[1, 1].set_xlabel("tSec (seconds)")
axs[1, 1].set_ylabel("range_rate (km/s)")

# # Show the plot
# plt.show()

fig.savefig("plots/meas_truth.pdf")
