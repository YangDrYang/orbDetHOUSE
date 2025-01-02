import os
import xml.etree.ElementTree as ET
import pandas as pd
from datetime import datetime, timedelta

norad_id = 46984  # NORAD ID
# Directory path where the .EOF files are located
directory = "refdata/"

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
eof_df.to_csv(directory + "od_ecef_id_" + str(norad_id) + ".csv", index=False)
