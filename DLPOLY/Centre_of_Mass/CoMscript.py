#This script will extract your molecule, calculate the centre of mass for it,
#and create a new file with the trajectory plot for your whole molecule + a new CoM atom.
#you can extract the CoM atoms only if you want (called 'X'). Note the final config needs a header and first config lines

import numpy as np

def extract_5FU_molecules(input_file, temp_file, expected_sequence, num_atoms):
    """
    Extracts 5FU molecules from the input file and writes them to a temporary file.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(temp_file, 'w') as outfile:
        for i in range(len(lines) - 2 * num_atoms + 1):
            sequence = [
                lines[i + 2 * j][0]  # First character of the line
                for j in range(num_atoms)
            ]
            if sequence == expected_sequence:
                outfile.writelines(lines[i:i + 2 * num_atoms])

def read_atoms_from_file(filename):
    atoms = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):
            atom_info = lines[i].split()
            coordinates = lines[i + 1].split()
            atom_type = atom_info[0]
            mass = float(atom_info[1])
            charge = float(atom_info[2])
            x, y, z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
            atoms.append([atom_type, mass, charge, x, y, z])
    return atoms

def calculate_center_of_mass(atoms):
    masses = np.array([atom[1] for atom in atoms])
    coordinates = np.array([[atom[3], atom[4], atom[5]] for atom in atoms])
    total_mass = np.sum(masses)
    weighted_coordinates = np.dot(masses, coordinates)
    center_of_mass = weighted_coordinates / total_mass
    return center_of_mass

def write_atoms_to_file(filename, atoms):
    with open(filename, 'w') as file:
        for atom in atoms:
            file.write(f"{atom[0]} {atom[1]} {atom[2]}\n")
            file.write(f"{atom[3]} {atom[4]} {atom[5]}\n")

def process_molecule_group(atoms, molecule_size):
    molecule_centers = []
    updated_atoms = []
    for i in range(0, len(atoms), molecule_size):
        molecule_atoms = atoms[i:i + molecule_size]
        center_of_mass = calculate_center_of_mass(molecule_atoms)
        dummy_atom = ['X', 0, 0, *center_of_mass]
        molecule_atoms.append(dummy_atom)
        updated_atoms.extend(molecule_atoms)
        molecule_centers.append(center_of_mass)
    return molecule_centers, updated_atoms

def main():
    # Input and output file paths
    input_file = 'history_input.txt'
    temp_file = 'temp_output.txt'
    output_file = 'history_output.txt'
    
    # Parameters for 5FU molecule
    expected_sequence = ['N', 'C', 'N', 'C', 'C', 'C', 'H', 'O', 'H', 'O', 'H', 'F']
    num_atoms = len(expected_sequence)
    molecule_size = 12

    # Step 1: Extract 5FU molecules
    print("Extracting 5FU molecules...")
    extract_5FU_molecules(input_file, temp_file, expected_sequence, num_atoms)
    print(f"5FU molecules extracted to {temp_file}.")

    # Step 2: Read extracted atoms
    print("Reading extracted atoms...")
    atoms = read_atoms_from_file(temp_file)

    # Step 3: Process molecules and calculate centers of mass
    print("Processing molecules...")
    molecule_centers, updated_atoms = process_molecule_group(atoms, molecule_size)

    # Step 4: Write updated atoms to output file
    print("Writing updated atoms to output file...")
    write_atoms_to_file(output_file, updated_atoms)

    # Output results
    print(f"Centers of Mass for all molecules: {molecule_centers}")
    print(f"Updated atomic data has been written to {output_file}")

if __name__ == "__main__":
    main()
