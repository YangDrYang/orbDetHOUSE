# Generate initial orbit via TLE

from downloadTLE import tleDownload
from skyfield.api import load as skyfieldLoad
import yaml

import pandas as pd

# Read the CSV file into a pandas DataFrame
obsDF = pd.read_csv("ccdata/meas_data.csv")

pathName = "pyscripts/"
yamlName = "InputStarlink.yml"
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
rECIFirstMeasEpoch = rvFirstMeasEpoch.position.km
vECIFirstMeasEpoch = rvFirstMeasEpoch.velocity.km_per_s

print(rECIFirstMeasEpoch)
print(vECIFirstMeasEpoch)
