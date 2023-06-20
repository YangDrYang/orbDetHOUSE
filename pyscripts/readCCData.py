import os
import pandas as pd
import datetime as dt

# Directory path
directory_path = "ccdata/20210522+23_LEO/"
# File name
file_name = "MatchedPasses600_20210522+23.Dat"
# Define the file path
file_path = directory_path + file_name

# Read the file into a DataFrame
df = pd.read_csv(file_path, delimiter="\s+", header=None)

# Filter rows based on the condition
filtered_df = df[df[3] == 45767]

# Retrieve the corresponding first column
file_directories = filtered_df[0]

# Print the file directories
print(file_directories)

# Extract file names from the file directories
file_names = [os.path.basename(path.replace("\\", "/")) for path in file_directories]

# Print the file names
print(file_names)


# Create an empty DataFrame to store the values
data_df = pd.DataFrame()

# Iterate over the file names
for file_name in file_names:
    file_path = os.path.join(directory_path, file_name)
    # Read the file and append its contents to the result DataFrame
    data = pd.read_csv(file_path, header=None, delimiter="\s+")
    data_df = pd.concat([data_df, data], ignore_index=True)

# Print the result DataFrame
print(data_df)

# Create a new DataFrame to store the modified data
obs_df = pd.DataFrame()

# Iterate over the rows of the original DataFrame
for index, row in data_df.iterrows():
    # Extract the date and time columns
    year, month, day, hour, minute, second = (
        int(row[0]),
        int(row[1]),
        int(row[2]),
        int(row[3]),
        int(row[4]),
        row[5],
    )

    # Convert the date and time to datetime object
    date_time = dt.datetime(year, month, day, hour, minute)

    # Calculate the Modified Julian Date (MJD)
    mjd = ((date_time - dt.datetime(1858, 11, 17)).total_seconds() + second) / (
        24 * 60 * 60
    )

    # Extract the angle columns
    ra, dec = row[6], row[7]

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
print(obs_df)

obs_df.to_csv("ccdata/meas_data.csv", index=False)
