import re

"""
...
IMPORTANT:
- Do NOT overly rely on this script for water removal.
- Best practice is to match simulation box dimensions with ice seed sizes 
  in two directions to ensure proper periodicity.

Script functionality:
Process a Gromacs .gro file by retaining the ice region and removing water molecules 
located inside the ice_zone or outside the simulation box.
When using gmx solvate to add water molecules, some waters may appear inside the ice region,
which should be removed. Also, some waters outside the box can be removed to tidy the system.

Parameters:
- input_file: input .gro file path/name
- output_file: output .gro file path/name
- box_dimensions: manually set box size [x, y, z] (units: nm)
- z_buffer: buffer size in z direction to adjust ice_zone boundaries

Note:
After deleting waters, atom numbering might become non-continuous, but this won't affect 
subsequent simulations. It is recommended to run an energy minimization (EM) step after 
and use the output file for further simulations.
"""

def process_gro_file(input_file, output_file, box_dimensions, z_buffer=0.0):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Header lines
    header = lines[:2]
    box_line = lines[-1]
    atoms = lines[2:-1]

    # Parse atoms and separate into fSOL and SOL
    fsol_atoms = []
    sol_atoms = []
    for line in atoms:
        residue_name = line[5:10].strip()
        if residue_name == "fSOL":
            fsol_atoms.append(line)
        elif residue_name == "SOL":
            sol_atoms.append(line)
    
    # Get z-coordinate range for fSOL
    z_coords = [float(line[36:44].strip()) for line in fsol_atoms]
    z_min, z_max = min(z_coords), max(z_coords)

    # Define the ice zone with user-defined box dimensions
    ice_zone = {
        "x_min": 0.0,
        "x_max": box_dimensions[0],
        "y_min": 0.0,
        "y_max": box_dimensions[1],
        "z_min": max(0.0, z_min - z_buffer),  # z_min adjusted by buffer
        "z_max": min(box_dimensions[2], z_max + z_buffer),  # z_max adjusted by buffer
    }

    # Check and filter SOL atoms
    filtered_sol_atoms = []
    deleted_molecule_count = 0
    molecule_start = None
    for i, line in enumerate(sol_atoms):
        if i % 4 == 0:  # Start of a new water molecule
            molecule_start = i
        
        # Extract coordinates
        x, y, z = map(float, (line[20:28].strip(), line[28:36].strip(), line[36:44].strip()))
        
        # Check if atom is in ice_zone or out of box dimensions
        in_ice_zone = (
            ice_zone["x_min"] <= x <= ice_zone["x_max"] and
            ice_zone["y_min"] <= y <= ice_zone["y_max"] and
            ice_zone["z_min"] <= z <= ice_zone["z_max"]
        )
        out_of_box = (
            x < 0.0 or x > box_dimensions[0] or
            y < 0.0 or y > box_dimensions[1] or
            z < 0.0 or z > box_dimensions[2]
        )

        # If any atom of the molecule is inside ice_zone or outside box, skip whole molecule
        if in_ice_zone or out_of_box:
            if i % 4 == 3:  # End of molecule
                deleted_molecule_count += 1
            continue

        # If molecule is valid, keep all its atoms
        if i % 4 == 3:  # End of water molecule
            filtered_sol_atoms.extend(sol_atoms[molecule_start:molecule_start + 4])
    
    # Renumber atoms
    combined_atoms = fsol_atoms + filtered_sol_atoms
    renumbered_atoms = []
    for atom_index, line in enumerate(combined_atoms, start=1):
        new_line = line[:15] + f"{atom_index:>5}" + line[20:]
        renumbered_atoms.append(new_line)

    # Update total number of atoms
    total_atoms = len(renumbered_atoms)
    header[1] = f"{total_atoms}\n"

    # Write output file
    with open(output_file, 'w') as file:
        file.writelines(header)
        file.writelines(renumbered_atoms)
        file.write(box_line)
    
    print(f"Processing complete! Removed {deleted_molecule_count} water molecules. Output file: {output_file}")


# ============= Configuration =============
# Input/output file paths - suggest placing script and input in same folder
input_gro = "ice-water.gro"          # Input file
output_gro = "ice-water-output.gro"  # Output file

# Manually set box dimensions (units: nm)
box_dimensions = [4.0, 4.0, 10.0]  # x, y, z dimensions in gro file

# z direction buffer size (units: nm)
z_buffer = 0.1  # adjust as needed, default 0.1

# Run main function
process_gro_file(input_gro, output_gro, box_dimensions, z_buffer)
