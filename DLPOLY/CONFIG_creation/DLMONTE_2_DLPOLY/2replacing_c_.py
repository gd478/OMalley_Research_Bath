input_file = "dlmonte.CONFIG"
output_file = "dlpoly.CONFIG"

# Read the content of the input file
with open(input_file, "r") as file:
    content = file.read()

# Replace occurrences of 'c' with incrementing numbers
counter = 1
def replace_c(match):
    global counter
    replacement = str(counter)
    counter += 1
    return replacement

import re
new_content = re.sub(r'c', replace_c, content)

# Write the modified content to the output file
with open(output_file, "w") as file:
    file.write(new_content)

print(f"Processed file saved as {output_file}")
