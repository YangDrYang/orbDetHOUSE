import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
from scipy.stats import skew, kurtosis
import seaborn as sns

# Directory path
out_folder_path = "out/out_ccdata/"
plot_folder_path = "plots/"

norad_id = 46984
meas_file = "ccdata/meas_data_id_" + str(norad_id) + ".csv"
stn_file = "ccdata/stn_eci_coordinates.csv"
# od_ref_data_file = "refdata/od_ref_id_" + str(norad_id) + ".csv"
od_ref_data_file = "refdata/od_ref_id_" + str(norad_id) + "_from_cpf.csv"

# # *************** pre-residuals


# pre measurement residuals
pre_res_file = "ccdata/ccdata_pre_res_id_" + str(norad_id) + ".csv"
pre_res_plot_file = (
    plot_folder_path + "ccdata_pre_res_plot_id_" + str(norad_id) + ".pdf"
)
pre_res_density_plot_file = (
    plot_folder_path + "ccdata_pre_res_density_id_" + str(norad_id) + ".pdf"
)
pre_res_all_plot_file = (
    plot_folder_path + "ccdata_pre_res_all_plot_id_" + str(norad_id) + ".pdf"
)


processing.process_pre_res_ccdata(
    od_ref_data_file, stn_file, meas_file, pre_res_file, []
)

pre_res_df = pd.read_csv(pre_res_file)
print(pre_res_df)
dates = pd.to_datetime(
    pre_res_df.loc[150:, "MJD"] + 2400000.5, unit="D", origin="julian"
)

ra_residuals = pre_res_df["RA"]  # Assuming this is your RA residuals column

ra_skew = skew(pre_res_df["RA"])
ra_kurt = kurtosis(pre_res_df["RA"])
dec_skew = skew(pre_res_df["Dec"])
dec_kurt = kurtosis(pre_res_df["Dec"])

<<<<<<< HEAD
# Set up the figure and subplots with a 1x2 grid
fig, axs = plt.subplots(1, 2, figsize=(6, 3))

# Plot histograms for ra_residuals and dec_residuals
sns.histplot(
    ra_residuals,
    bins=50,
    kde=True,
    color="blue",
    label="RA Residuals (arcsec)",
    ax=axs[0],
)
sns.histplot(
    dec_residuals,
    bins=50,
    kde=True,
    color="green",
    label="Dec Residuals (arcsec)",
    ax=axs[1],
)

# Set titles, labels, and legends
axs[0].set_title("RA Residuals")
axs[0].set_ylabel("Density")
axs[0].legend()
axs[1].set_title("Dec Residuals")
axs[1].set_ylabel("Density")
axs[1].legend()

# Adjust the spacing between subplots
plt.tight_layout()
# Save the figure
plt.savefig(pre_res_density_plot_file)

ra_mean = np.mean(ra_residuals)
ra_std = np.std(ra_residuals)
dec_mean = np.mean(dec_residuals)
dec_std = np.std(dec_residuals)
ra_rms = np.sqrt(np.mean(np.square(ra_residuals)))
dec_rms = np.sqrt(np.mean(np.square(dec_residuals)))
ra_skew = skew(ra_residuals)
ra_kurt = kurtosis(ra_residuals)
dec_skew = skew(dec_residuals)
dec_kurt = kurtosis(dec_residuals)

print("right ascension mean (in arcseconds):    ", ra_mean)
print("right ascension std (in arcseconds):    ", ra_std)
print("right ascension rms (in arcseconds):    ", ra_rms)
=======
>>>>>>> parent of 8b55a691 (corrected cc stn coor and updated pre-residuals)
print("right ascension skewness:    ", ra_skew)
print("right ascension kurtosis:    ", ra_kurt)
print("declination skewness:    ", dec_skew)
print("declination kurtosis:    ", dec_kurt)


# Create a new figure and axes for plot
fig, ax1 = plt.subplots(figsize=(6, 8))
# Set the font size
plt.rcParams.update({"font.size": 12})
# Plot RA vs MJD
ax2 = ax1.twinx()

ax1.scatter(
    dates,
    np.degrees(pre_res_df.loc[150:, "RA"]) * 3600,
    color="blue",
    s=5,
    label="RA Residuals",
)
ax1.set_ylabel("RA Residuals \n(arcsec)", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")
# ax1.set_ylim(-5e-6, 7.5e-6)
# ax1.set_ylim(-30, 40)

# add 20 mins shift to distinguish ra and dec
dates = pd.to_datetime(
    pre_res_df.loc[150:, "MJD"] + 2400000.5 + 20 / 1440, unit="D", origin="julian"
)
ax2.scatter(
    dates,
<<<<<<< HEAD
    dec_residuals,
    color="green",
=======
    np.degrees(pre_res_df.loc[150:, "Dec"]) * 3600,
    color="red",
>>>>>>> parent of 8b55a691 (corrected cc stn coor and updated pre-residuals)
    s=5,
    label="Dec Residuals",
)
ax2.set_ylabel("Dec Residuals \n(arcsec)", color="green")
ax2.tick_params(axis="y", labelcolor="green")
# ax2.set_ylim(-25, 40)


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

# Rotate x-axis tick labels and adjust label spacing
fig.autofmt_xdate(rotation=45)
plt.tight_layout()

# Save the figure
plt.savefig(pre_res_plot_file)


# Generate only one figure

# Create a new figure with a 2x1 grid
fig = plt.figure(figsize=(8, 6))

# Subplot 1: Scatter plot for RA and Dec Residuals (Top)
ax1 = fig.add_subplot(2, 1, 1)
ax2 = ax1.twinx()

ax1.scatter(
    dates,
    ra_residuals,
    color="blue",
    s=5,
    label="RA Residuals",
)
ax1.set_ylabel("RA Residuals \n(arcsec)", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")

# add 30 mins shift to distinguish ra and dec
dates = pd.to_datetime(
    pre_res_df.loc[:, "MJD"] + 2400000.5 + 30 / 1440, unit="D", origin="julian"
)
ax2.scatter(
    dates,
    dec_residuals,
    color="green",
    s=5,
    label="Dec Residuals",
)
ax2.set_ylabel("Dec Residuals \n(arcsec)", color="green")
ax2.tick_params(axis="y", labelcolor="green")

# # Format x-axis as hours
# hours = mdates.HourLocator(interval=6)
# hour_format = mdates.DateFormatter("%H:%M")
# ax1.xaxis.set_major_locator(hours)
# ax1.xaxis.set_major_formatter(hour_format)

# # Show only the date for the first epoch of each day
# dates_first_epoch = np.unique([date.date() for date in dates])
# for date in dates_first_epoch:
#     ax1.axvline(date, color="gray", linestyle="--", alpha=0.5)
#     ax2.axvline(date, color="gray", linestyle="--", alpha=0.5)

# # Show only 1 hour before and 1 hour after the data on x-axis
# start_time = min(dates) - dt.timedelta(hours=1)
# end_time = max(dates) + dt.timedelta(hours=1)
# ax1.set_xlim(start_time, end_time)

# # Rotate x-axis tick labels and adjust label spacing
# fig.autofmt_xdate(rotation=45)

start_time = min(dates) - dt.timedelta(hours=1)
end_time = max(dates) + dt.timedelta(hours=2)

# Show only the date for the first epoch of each day
dates_first_epoch = np.unique([date.date() for date in dates])
for date in dates_first_epoch[1:]:
    ax1.axvline(date, color="gray", linestyle="--", alpha=0.5)
    ax2.axvline(date, color="gray", linestyle="--", alpha=0.5)

# Set x-axis ticks and labels for the desired time range
time_ticks = pd.date_range(
    start=start_time, end=end_time, freq="180T"
)  # Adjust frequency as needed
time_tick_labels = [time.strftime("%H:%M") for time in time_ticks]

# Set x-axis ticks and labels for ax1 and ax2
ax1.set_xticks(time_ticks)
ax1.set_xticklabels(time_tick_labels, rotation=45, ha="right")
ax2.set_xticks(time_ticks)
ax2.set_xticklabels(time_tick_labels, rotation=45, ha="right")

# Set the title for the combined plot
fig.suptitle("RA and Dec Residuals", fontsize=16)

# # Rotate x-axis tick labels and adjust label spacing in the scatter plot
# fig.autofmt_xdate(rotation=45)
# ax1.set_title("RA and Dec Residuals")

# Subplots 2 and 3: Density histograms for RA and Dec Residuals (Bottom)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

sns.histplot(
    ra_residuals,
    bins=50,
    kde=True,
    color="blue",
    label="RA Residuals",
    ax=ax3,
)

ax3.set_ylabel("Density")
ax3.set_title("RA Residuals Distribution")

sns.histplot(
    dec_residuals,
    bins=50,
    kde=True,
    color="green",
    label="Dec Residuals",
    ax=ax4,
)

ax4.set_ylabel("Density")
ax4.set_title("Dec Residuals Distribution")

# Adjust the spacing between subplots
plt.tight_layout()

# Save the figure
plt.savefig(pre_res_all_plot_file)
