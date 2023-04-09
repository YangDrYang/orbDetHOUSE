import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

# Read the CSV file into a Pandas DataFrame
df = pd.read_csv("out_/trajectory_truth.csv")

# Extract the x, y, and z columns from the DataFrame
x = df["x"]
y = df["y"]
z = df["z"]

# Create a 3D plot using Matplotlib and Seaborn
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")
ax.scatter(xs=x, ys=y, zs=z, s=50, c=z, cmap="viridis", edgecolor="none")

# Add labels and a title to the plot
ax.set_xlabel("X Label", fontsize=14)
ax.set_ylabel("Y Label", fontsize=14)
ax.set_zlabel("Z Label", fontsize=14)
ax.set_title("3D Scatter Plot", fontsize=16)

# Customize the legend
legend = ax.legend(title="Z Values", fontsize=12)
legend.get_title().set_fontsize(12)

# Show the plot
plt.show()
