input_file = "dlmonte.CONFIG"
output_file = "dlpoly.CONFIG"

# Read the content of the input file
with open(input_file, "r") as file:
    lines = file.readlines()

# Remove lines containing 'MOLECULE' or 'NUMMOL'
filtered_lines = [line for line in lines if 'MOLECULE' not in line and 'NUMMOL' not in line]

# Write the modified content to the output file
with open(output_file, "w") as file:
    file.writelines(filtered_lines)

print(f"Processed file saved as {output_file}")