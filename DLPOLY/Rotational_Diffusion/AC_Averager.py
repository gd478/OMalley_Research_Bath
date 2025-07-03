#gd478 AC h5 file averager#
import os
import glob
import h5py
import numpy as np
import pandas as pd

# Get the current working directory
directory = os.getcwd()

# Get a list of all .h5 files in the current directory
h5_files = glob.glob(os.path.join(directory, "*.h5"))

# Initialize lists to accumulate time and AC data across all files
time_data = None
ac_data_list = []

# Iterate over each .h5 file
for h5_file in h5_files:
    with h5py.File(h5_file, 'r') as file:
        dataset_names = list(file.keys())
        if len(dataset_names) == 0:
            raise ValueError(f"No datasets found in the HDF5 file {h5_file}.")

        # Assume the first dataset contains the desired AC data
        dataset_name = dataset_names[0]
        data = file[dataset_name][:]

        # Check if the "time" variable is present
        variable_name = "time"
        if variable_name in file:
            time = file[variable_name][:]
            if time_data is None:
                # Initialize time_data only once from the first file
                time_data = time
            elif not np.array_equal(time_data, time):
                raise ValueError("Time arrays are not identical across all files.")
        else:
            raise ValueError(f"'time' dataset is missing in file {h5_file}.")

        # Accumulate AC data for averaging
        ac_data_list.append(data)

# Convert list of AC data arrays to a 2D NumPy array and calculate the mean across files
ac_data_stack = np.stack(ac_data_list, axis=0)
ac_data_avg = np.mean(ac_data_stack, axis=0)

# Combine time and averaged AC data into a DataFrame for CSV export
combined_data = np.column_stack((time_data, ac_data_avg))
df = pd.DataFrame(combined_data, columns=['time (ps)', 'AC avg'])

# Define the output CSV file path
output_csv = os.path.join(directory, 'final.csv')

# Save the DataFrame to CSV
df.to_csv(output_csv, index=False)
print(f"Successfully saved the averaged data to {output_csv}")
