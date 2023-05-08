import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import processing
import seaborn as sns
import numpy as np

filters = ["house", "ukf", "cut4", "cut6"]

# folder_path = "out_/"
folder_path = "out/"

# Create a new figure and axes for logarithmic plot
fig, ax = plt.subplots(ncols=2, figsize=(10, 5))

run_times = []
for filter_type in filters:
    # Read CSV file into a pandas dataframe
    df = pd.read_csv(folder_path + "run_times_" + filter_type + ".csv")

    run_times.append(df[filter_type])
