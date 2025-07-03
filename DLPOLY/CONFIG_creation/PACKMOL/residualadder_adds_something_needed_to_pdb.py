input_file = "water.pdb"   # Replace with the path to your input PDB file
output_file = "waterr.pdb"  # Output file to save the modified content

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if len(line) >= 24:
            # Replace character at position 24 with '1' without shifting other characters
            modified_line = line[:23] + "1" + line[24:]
        else:
            # Leave the line unchanged if it's shorter than 24 characters
            modified_line = line
        outfile.write(modified_line)

print(f"Modified PDB file has been saved as {output_file}")
