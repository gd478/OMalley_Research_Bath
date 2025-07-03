#arrange your data in an excel file called test.xlsx
#column 1 = the x values of your data. Column 2 = the y values
#column 3 = the desired x values that you would like to interpolate new y values for
#run script: python3 interpolater.py test.xlsx

import pandas as pd
from scipy import interpolate

# Read the Excel file
file_path = 'test.xlsx'
df = pd.read_excel(file_path)

# Extract the columns
x_values = df.iloc[:, 0].values  # First column (x values for original data)
y_values = df.iloc[:, 1].values  # Second column (y values for original data)
target_x_values = df.iloc[:, 2].values  # Third column (target x values)

# Create a cubic interpolation function (better for non-linear data)
# You can also use 'PchipInterpolator' or 'CubicSpline' if preferred
interp_function = interpolate.interp1d(x_values, y_values, kind='cubic', fill_value="extrapolate")

# Use the interpolation function to generate y values for the target x values
interpolated_y_values = interp_function(target_x_values)

# Add the interpolated y values to the dataframe as a new column (4th column)
df['Interpolated Y'] = interpolated_y_values

# Save the new dataframe with the interpolated values to a new Excel file
output_file_path = 'test_with_interpolated_y.xlsx'
df.to_excel(output_file_path, index=False)

print(f"Interpolation complete. The new file is saved as '{output_file_path}'.")
