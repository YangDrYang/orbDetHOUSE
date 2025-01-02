import pandas as pd
import matplotlib.pyplot as plt

# Load prop_results.csv into a dataframe
prop_results_path = "out/out_prop/prop_results.csv"
prop_results_df = pd.read_csv(prop_results_path)

# Load od_eci_id_46984.csv into a dataframe
od_eci_path = "refdata/od_eci_id_46984.csv"
od_eci_df = pd.read_csv(od_eci_path)
od_eci_df = od_eci_df.iloc[0:1561, 0:7]
od_eci_df.columns = ["mjd", "x", "y", "z", "vx", "vy", "vz"]
print(od_eci_df)


# Extract every 6th row from prop_results_df
prop_results_extracted_df = prop_results_df.iloc[::6]

# Reset the index of prop_results_extracted_df to match od_eci_df
prop_results_extracted_df.reset_index(drop=True, inplace=True)

# Drop the first column from od_eci_df and prop_results_extracted_df
od_eci_extracted_df = od_eci_df.iloc[:, 1:7]
prop_results_extracted_df = prop_results_extracted_df.iloc[:, 1:7]
print(prop_results_extracted_df)

# Rename the columns in od_eci_df
od_eci_extracted_df.columns = ["x", "y", "z", "vx", "vy", "vz"]
print(od_eci_extracted_df)


# Subtract prop_results_extracted_df from od_eci_df
differences_df = od_eci_extracted_df.subtract(prop_results_extracted_df)

# Add the first column of differences_df equal to MJD
differences_df.insert(0, "mjd", od_eci_df["mjd"])
print(differences_df)

# Plot x, y, z, vx, vy, vz versus MJD in separate subfigures
fig, axs = plt.subplots(nrows=3, figsize=(6, 8))

# Plot x, y, z, vx, vy, vz versus MJD in separate subplots
components = ["x", "y", "z"]
for i, component in enumerate(components):
    axs[i].plot(differences_df["mjd"], differences_df[component])
    axs[i].set_xlabel("mjd")
    axs[i].set_ylabel(f"difference in {component}")
    axs[i].set_title(f"differences in {component} versus mjd")

plt.tight_layout()
plt.show()

# Plot x, y, z, vx, vy, vz versus MJD in separate subfigures
fig, axs = plt.subplots(nrows=3, figsize=(6, 8))

# Plot x, y, z, vx, vy, vz versus MJD in separate subplots
components = ["vx", "vy", "vz"]
for i, component in enumerate(components):
    axs[i].plot(differences_df["mjd"], differences_df[component])
    axs[i].set_xlabel("mjd")
    axs[i].set_ylabel(f"difference in {component}")
    axs[i].set_title(f"differences in {component} versus mjd")

plt.tight_layout()
plt.show()
