import pandas as pd
import matplotlib.pyplot as plt
import processing
import seaborn as sns


# filters = ["ukf"]
# filters = ["cut6"]
filters = ["house", "ukf", "cut4", "cut6"]
# filters = ["house", "ukf"]

folder_path = "out_/"


for filter_type in filters:
    processing.process_err_each_filter(filter_type, folder_path)
