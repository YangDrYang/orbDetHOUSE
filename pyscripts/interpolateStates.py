# Define a namedtuple to hold the x and y values
from collections import namedtuple
import pandas as pd
import sys

norad_id = 46984
directory_path = "ccdata/"  # Directory path
file_name = "meas_data_id_" + str(norad_id) + ".csv"  # File name
file_path = directory_path + file_name
# Read the measurement CSV file into a pandas DataFrame
meas_df = pd.read_csv(file_path)

# Retrieving the command-line argument
if len(sys.argv) > 1:
    orbECIFile = sys.argv[1]
else:
    orbECIFile = input("Enter the orbit eci file: ")
# Read the od reference (in ECI) CSV file into a pandas DataFrame
# od_eci_df = pd.read_csv("refdata/od_eci_id_" + str(norad_id) + ".csv")
od_eci_df = pd.read_csv(orbECIFile)
# Check if the last column contains the keyword "Unnamed"
last_column_name = od_eci_df.columns[-1]
if "Unnamed" in last_column_name:
    # Drop the last column from the DataFrame
    od_eci_df = od_eci_df.iloc[:, :-1]
print(od_eci_df)


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


file_without_extension = orbECIFile.split(".")[0]
last_three_chars = file_without_extension[-3:]
if last_three_chars != "cpf":
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

    # Iterate over each row in the measurement dataframe
    for index, row in meas_df.iterrows():
        # Extract the target MJD for interpolation
        target_mjd = row["MJD"]

        # Extract the nearest 10 points for interpolation
        nearest_points = od_eci_df.iloc[
            (od_eci_df["MJD"] - target_mjd).abs().argsort()[:nOrd]
        ]

        # Create DataPoint objects for the nearest points
        data_points_x = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["X_ECI"])
        ]
        data_points_y = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["Y_ECI"])
        ]
        data_points_z = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["Z_ECI"])
        ]
        data_points_vx = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["VX_ECI"])
        ]
        data_points_vy = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["VY_ECI"])
        ]
        data_points_vz = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["VZ_ECI"])
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
    od_ref_df.to_csv("refdata/od_ref_id_" + str(norad_id) + ".csv", index=False)

else:
    od_ref_df = pd.DataFrame(
        columns=[
            "MJD",
            "Interpolated_X",
            "Interpolated_Y",
            "Interpolated_Z",
        ]
    )

    # Iterate over each row in the measurement dataframe
    for index, row in meas_df.iterrows():
        # Extract the target MJD for interpolation
        target_mjd = row["MJD"]

        # Extract the nearest 10 points for interpolation
        nearest_points = od_eci_df.iloc[
            (od_eci_df["MJD"] - target_mjd).abs().argsort()[:nOrd]
        ]

        # Create DataPoint objects for the nearest points
        data_points_x = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["X_ECI"])
        ]
        data_points_y = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["Y_ECI"])
        ]
        data_points_z = [
            DataPoint(x, y)
            for x, y in zip(nearest_points["MJD"], nearest_points["Z_ECI"])
        ]

        # Perform the Lagrange interpolation for X/Y/Z/VX/VY/VZ
        interpolated_x = interpolate(data_points_x, target_mjd, len(data_points_x))
        interpolated_y = interpolate(data_points_y, target_mjd, len(data_points_y))
        interpolated_z = interpolate(data_points_z, target_mjd, len(data_points_z))

        # Create a new row with the interpolated values
        new_row = {
            "MJD": target_mjd,
            "Interpolated_X": interpolated_x,
            "Interpolated_Y": interpolated_y,
            "Interpolated_Z": interpolated_z,
        }

        # Concatenate the new row to the interpolated dataframe
        od_ref_df = pd.concat([od_ref_df, pd.DataFrame([new_row])], ignore_index=True)

    # # Print the interpolated dataframe
    # print(od_ref_df)

    # Save od_ref_df to a CSV file
    od_ref_df.to_csv(
        "refdata/od_ref_id_" + str(norad_id) + "_from_cpf.csv", index=False
    )
