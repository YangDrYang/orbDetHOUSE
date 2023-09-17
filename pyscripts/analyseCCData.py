import pandas as pd
import matplotlib.pyplot as plt
import processing
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
from scipy.stats import skew, kurtosis

# Directory path
out_folder_path = "out/out_ccdata/"
plot_folder_path = "plots/"

norad_id = 46984
meas_file = "ccdata/meas_data_id_" + str(norad_id) + ".csv"
stn_file = "ccdata/stn_eci_coordinates.csv"
od_ref_data_file = "refdata/od_ref_id_" + str(norad_id) + ".csv"

# # *************** pre-residuals


# pre measurement residuals
pre_res_file = "ccdata/ccdata_pre_res_id_" + str(norad_id) + ".csv"
pre_res_plot_file = plot_folder_path + "ccdata_pre_res_id_" + str(norad_id) + ".pdf"

processing.process_pre_res_each_filter_ccdata(
    od_ref_data_file, stn_file, meas_file, pre_res_file, []
)

pre_res_df = pd.read_csv(pre_res_file)
dates = pd.to_datetime(pre_res_df["MJD"] + 2400000.5, unit="D", origin="julian")

ra_residuals = pre_res_df["RA"]  # Assuming this is your RA residuals column

ra_skew = skew(pre_res_df["RA"])
ra_kurt = kurtosis(pre_res_df["RA"])
dec_skew = skew(pre_res_df["Dec"])
dec_kurt = kurtosis(pre_res_df["Dec"])

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
    np.degrees(pre_res_df["RA"]) * 3600,
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
    pre_res_df["MJD"] + 2400000.5 + 20 / 1440, unit="D", origin="julian"
)
ax2.scatter(
    dates,
    np.degrees(pre_res_df["Dec"]) * 3600,
    color="red",
    s=5,
    label="Dec Residuals",
)
ax2.set_ylabel("Dec Residuals \n(arcsec)", color="red")
ax2.tick_params(axis="y", labelcolor="red")
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
