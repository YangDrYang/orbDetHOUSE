import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

# Read the CSV file into a Pandas DataFrame
folder_path = "out_/"  # replace with the path to your folder
df = pd.read_csv(folder_path + "trajectory_truth.csv")

# Extract the x, y, and z columns from the DataFrame
x = df["x"] / 1000
y = df["y"] / 1000
z = df["z"] / 1000

# Create a 3D plot using Matplotlib and Seaborn
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")
ax.scatter(xs=x, ys=y, zs=z, s=50, c=z, cmap="viridis", edgecolor="none", linewidth=0.5)

# Add labels and a title to the plot
ax.set_xlabel("x component [km]", fontsize=14)
ax.set_ylabel("y component [km]", fontsize=14)
ax.set_zlabel("z component [km]", fontsize=14)
# ax.set_title("3D orbit scatter plot", fontsize=16)
plt.tight_layout()

# Customize the legend
# legend = ax.legend(title="Z Values", fontsize=12)
# legend.get_title().set_fontsize(12)

# # Show the plot
# plt.show()

fig.savefig("plots/traj_3d_truth.pdf")
