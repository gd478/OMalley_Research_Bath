import re

def replace_o_numbers(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    with open(output_file, 'w') as f:
        for line in lines:
            modified_line = re.sub(r'^(O|Si)\s+\d+', r'\1                c', line)
            modified_line = re.sub(r'(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)', r'\1 0', modified_line)
            f.write(modified_line)

# Example usage
input_filename = 'output.txt'  # Change this to your actual filename
output_filename = 'output2.txt'  # Output file
replace_o_numbers(input_filename, output_filename)
