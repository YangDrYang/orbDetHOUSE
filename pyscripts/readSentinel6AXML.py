import os
import xml.etree.ElementTree as ET
import pandas as pd
from datetime import datetime, timedelta

norad_id = 46984  # NORAD ID
# Directory path where the .EOF files are located
directory = "pyscripts/"

# Create an empty list to store data from all files
data = []

# Iterate over the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".EOF"):
        file_path = os.path.join(directory, filename)

        # Load the XML file
        tree = ET.parse(file_path)
        root = tree.getroot()

        # Find all OSV elements
        osv_elements = root.findall(".//OSV")

        # Iterate over each OSV element
        for osv in osv_elements:
            utc_str = osv.find(".//UTC").text.split("=")[1]
            utc = datetime.strptime(utc_str, "%Y-%m-%dT%H:%M:%S.%f")
            mjd = (utc - datetime(1858, 11, 17)) / timedelta(days=1)
            x = float(osv.find(".//X").text)
            y = float(osv.find(".//Y").text)
            z = float(osv.find(".//Z").text)
            vx = float(osv.find(".//VX").text)
            vy = float(osv.find(".//VY").text)
            vz = float(osv.find(".//VZ").text)

            # Add the values to the data list
            data.append([utc_str, mjd, x, y, z, vx, vy, vz])

# Create a dataframe from the data list
eof_df = pd.DataFrame(data, columns=["UTC", "MJD", "X", "Y", "Z", "VX", "VY", "VZ"])

# # Print the dataframe
# print(eof_df)

# Save eof_df to a CSV file
eof_df.to_csv(directory + "eof_data_" + str(norad_id) + ".csv", index=False)

directory_path = "ccdata/"  # Directory path
file_name = "meas_data_id_" + str(norad_id) + ".csv"  # File name
file_path = directory_path + file_name
# Read the measurement CSV file into a pandas DataFrame
meas_df = pd.read_csv(file_path)

# Create an empty dataframe to store the interpolated values
od_ref_df = pd.DataFrame(
    columns=[
        "MJD",
        "Interpolated_X",
        "Interpolated_Y",
        "Interpolated_Z",
        "Interpolated_VX",
        "Interpolated_VY",
        "Interpolated_VZ",
    ]
)


# Define a namedtuple to hold the x and y values
from collections import namedtuple

DataPoint = namedtuple("DataPoint", ["x", "y"])

# Order number
nOrd = 10


# Function to interpolate using Lagrange's formula
def interpolate(f, xi, n):
    result = 0.0
    for i in range(n):
        term = f[i].y
        for j in range(n):
            if j != i:
                term *= (xi - f[j].x) / (f[i].x - f[j].x)
        result += term
    return result


# Iterate over each row in the measurement dataframe
for index, row in meas_df.iterrows():
    # Extract the target MJD for interpolation
    target_mjd = row["MJD"]

    # Extract the nearest 10 points for interpolation
    nearest_points = eof_df.iloc[(eof_df["MJD"] - target_mjd).abs().argsort()[:nOrd]]

    # Create DataPoint objects for the nearest points
    data_points_x = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["X"])
    ]
    data_points_y = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["Y"])
    ]
    data_points_z = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["Z"])
    ]
    data_points_vx = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["VX"])
    ]
    data_points_vy = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["VY"])
    ]
    data_points_vz = [
        DataPoint(x, y) for x, y in zip(nearest_points["MJD"], nearest_points["VZ"])
    ]

    # Perform the Lagrange interpolation for X/Y/Z/VX/VY/VZ
    interpolated_x = interpolate(data_points_x, target_mjd, len(data_points_x))
    interpolated_y = interpolate(data_points_y, target_mjd, len(data_points_y))
    interpolated_z = interpolate(data_points_z, target_mjd, len(data_points_z))
    interpolated_vx = interpolate(data_points_vx, target_mjd, len(data_points_vx))
    interpolated_vy = interpolate(data_points_vy, target_mjd, len(data_points_vy))
    interpolated_vz = interpolate(data_points_vz, target_mjd, len(data_points_vz))

    # Create a new row with the interpolated values
    new_row = {
        "MJD": target_mjd,
        "Interpolated_X": interpolated_x,
        "Interpolated_Y": interpolated_y,
        "Interpolated_Z": interpolated_z,
        "Interpolated_VX": interpolated_vx,
        "Interpolated_VY": interpolated_vy,
        "Interpolated_VZ": interpolated_vz,
    }

    # Concatenate the new row to the interpolated dataframe
    od_ref_df = pd.concat([od_ref_df, pd.DataFrame([new_row])], ignore_index=True)

# # Print the interpolated dataframe
# print(od_ref_df)


# Save od_ref_df to a CSV file
od_ref_df.to_csv(directory_path + "od_ref_data_" + str(norad_id) + ".csv", index=False)
