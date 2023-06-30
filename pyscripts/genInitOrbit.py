# Generate initial orbit via TLE

from downloadTLE import tleDownload
from skyfield.api import load as skyfieldLoad
import yaml

import pandas as pd

norad_id = 46984  # NORAD ID
directory_path = "ccdata/"  # Directory path
file_name = "meas_data_id_" + str(norad_id) + ".csv"  # File name
file_path = directory_path + file_name
# Read the CSV file into a pandas DataFrame
obsDF = pd.read_csv(file_path)

pathName = "yamls/"
yamlName = "inputSentinel6A.yml"
with open(pathName + yamlName) as f:
    inputDict = yaml.safe_load(f)

satDict = inputDict["satDict"]
outFile = inputDict["dictFile"]["tleOutFile"]
satdata = tleDownload(satDict, outFile)

ts = skyfieldLoad.timescale()
tProp = ts.utc(
    year=obsDF.iloc[0, 0],
    month=obsDF.iloc[0, 1],
    day=obsDF.iloc[0, 2],
    hour=obsDF.iloc[0, 3],
    minute=obsDF.iloc[0, 4],
    second=obsDF.iloc[0, 5],
)

print(obsDF.iloc[0, 0:6])

# obtain position and velocity at the next epoch from TLE
rvFirstMeasEpoch = satdata[0].at(tProp)
rECIFirstMeasEpoch = rvFirstMeasEpoch.position.km * 1000
vECIFirstMeasEpoch = rvFirstMeasEpoch.velocity.km_per_s * 1000

print(tProp)
print(rECIFirstMeasEpoch)
print(vECIFirstMeasEpoch)
