import pandas as pd
import matplotlib.pyplot as plt
import processing
import seaborn as sns
import sys

if len(sys.argv) < 2:
    folder_path = "out/out_dense/"
    # folder_path = "out/out_sparse/"
    # folder_path = "out/out_sparse_pearson/"
    # folder_path = "out/"
else:
    folder_path = sys.argv[1]

# filters = ["ukf"]
# filters = ["cut6"]
# filters = ["ukf", "cut4", "cut6", "house"]
filters = ["ukf", "cut4", "cut6", "house", "srhouse"]
# filters = ["house", "ukf"]

for filter_type in filters:
    processing.process_err_each_filter(filter_type, folder_path)
