import csv


def read_csv(file_path):
    data = []
    with open(file_path, "r") as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        for row in reader:
            data.append(
                [float(x) if x else 0.0 for x in row]
            )  # Convert empty strings to 0.0
    return data


def compute_differences(data1, data2):
    differences = []
    for row1, row2 in zip(data1, data2):
        diff_row = [row1[0]]  # Keep the timestamp
        diff_row += [x1 - x2 for x1, x2 in zip(row1[1:], row2[1:])]
        differences.append(diff_row)
    return differences


def write_csv(file_path, data):
    with open(file_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["tSec", "dx", "dy", "dz", "dvx", "dvy", "dvz"])  # Write header
        writer.writerows(data)


# File paths
file1 = "out/out_prop/out_propprop_results_py.csv"
file2 = "out/out_propprop_results_py.csv"
output_file = "out/out_prop/differences.csv"

# Read data from CSV files
data1 = read_csv(file1)
data2 = read_csv(file2)

# Compute differences
differences = compute_differences(data1, data2)

# Write differences to a new CSV file
write_csv(output_file, differences)

print(f"Differences have been written to {output_file}")
