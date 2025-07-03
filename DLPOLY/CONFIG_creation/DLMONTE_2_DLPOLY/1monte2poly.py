import re

def convert_dlmonte_to_dlpoly(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    counter = 1  # Initialize counter
    
    with open(output_file, 'w') as f:
        for line in lines:
            # Replace 'c' with an incrementing number, ensuring correct spacing
            modified_line = re.sub(r'^(O|Si|N|C_|C|H_N|O_|H|F)\s+c', lambda m: f'{m.group(1)}        {counter}', line)
            if re.search(r'^(O|Si|N|C_|C|H_N|O_|H|F)\s+c', line):
                counter += 1  # Increment counter only if a replacement is made
            
            # Remove the extra '0' at the end of coordinate lines
            modified_line = re.sub(r'(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\s+0', r'\1', modified_line)
            f.write(modified_line)

# Example usage
input_filename = 'dlmonte.CONFIG'  # Change this to your actual filename
output_filename = 'dlpoly.CONFIG'  # Output file
convert_dlmonte_to_dlpoly(input_filename, output_filename)
