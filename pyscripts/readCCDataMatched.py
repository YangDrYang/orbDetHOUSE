import os
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

# NORAD ID
norad_id = 46984
# Directory path
directory_path = "ccdata/"
# File name
file_name = "MatchedPasses_NoradID_46984.PAS"
# Define the file path
file_path = directory_path + file_name

# Read the file into a DataFrame
data_df = pd.read_csv(file_path, delimiter="\s+", header=None)

# Create a new DataFrame to store the modified data
obs_df = pd.DataFrame()

# Iterate over the rows of the original DataFrame
for index, row in data_df.iterrows():
    # Extract the date and time columns
    year, month, day, hour, minute, second = (
        int(row[1]),
        int(row[2]),
        int(row[3]),
        int(row[4]),
        int(row[5]),
        row[6],
    )

    # Convert the date and time to datetime object
    date_time = dt.datetime(year, month, day, hour, minute)

    # Calculate the Modified Julian Date (MJD)
    mjd = ((date_time - dt.datetime(1858, 11, 17)).total_seconds() + second) / (
        24 * 60 * 60
    )

    # Extract the angle columns
    ra, dec = row[7], row[8]

    # Create a new DataFrame with the modified values
    new_row = pd.DataFrame(
        {
            "Year": [year],
            "Month": [month],
            "Day": [day],
            "Hour": [hour],
            "Minute": [minute],
            "Second": [second],
            "MJD": [mjd],
            "RA": [ra],
            "Dec": [dec],
        }
    )

    # Concatenate the new row with the existing DataFrame
    obs_df = pd.concat([obs_df, new_row], ignore_index=True)

# # Print the new DataFrame
# print(obs_df)

# Sort the dataframe according to MJD in an acsending order
obs_df = obs_df.sort_values(by="MJD")

# Keep only obs_df rows that are before 2021/05/19
obs_df = obs_df[obs_df["MJD"] < 59353]

print(obs_df)

obs_df.to_csv("ccdata/meas_data_id_" + str(norad_id) + ".csv", index=False)


# Convert MJD to datetime objects
dates = obs_df["MJD"].apply(
    lambda mjd: dt.datetime(1858, 11, 17) + dt.timedelta(days=mjd)
)
# # Convert MJD to datetime objects
# dates = mdates.num2date(obs_df["MJD"])
# print(dates)

# Create subplots with shared x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 10))

# Plot RA vs MJD
ax1.scatter(dates, np.radians(obs_df["RA"]), color="blue")
ax1.set_ylabel("ra (radians)")
# ax1.set_title("ra vs date")

# Plot Dec vs MJD
ax2.scatter(dates, np.radians(obs_df["Dec"]), color="red")
ax2.set_xlabel("date (utc)")
ax2.set_ylabel("dec (radians)")
# ax2.set_title("dec vs date")

# Format x-axis as hours
hours = mdates.HourLocator(interval=5)
hour_format = mdates.DateFormatter("%H:%M")
ax2.xaxis.set_major_locator(hours)
ax2.xaxis.set_major_formatter(hour_format)

# Show only the date for the first epoch of each day
dates_first_epoch = np.unique([date.date() for date in dates])
print(dates_first_epoch)
for date in dates_first_epoch:
    ax1.axvline(date, color="gray", linestyle="--", alpha=0.5)
    ax2.axvline(date, color="gray", linestyle="--", alpha=0.5)

# Show only 1 hour before and 1 hour after the data on x-axis
start_time = min(dates) - dt.timedelta(hours=2)
end_time = max(dates) + dt.timedelta(hours=2)
ax2.set_xlim(start_time, end_time)


# Rotate x-axis tick labels and adjust label spacing
fig.autofmt_xdate(rotation=50)
ax2.tick_params(axis="x", rotation=45, labelsize=10)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.4)

plt.tight_layout()

# Save the figure
plt.savefig("ccdata/meas_data_id_" + str(norad_id) + ".png")
