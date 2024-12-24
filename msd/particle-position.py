import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非图形化的Agg后端
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
配合msd-particle-position.py使用
这个脚本生成的坐标位置
是另一个脚本的输入
'''

# Read PDB file and return coordinates for all models
def read_pdb(pdb_file):
    models = []
    with open(pdb_file, 'r') as f:
        model = []
        for line in f:
            if line.startswith("MODEL"):
                if model:
                    models.append(model)
                model = []
            if line.startswith("ATOM") or line.startswith("HETATM"):
                model.append(line)
        if model:
            models.append(model)
    return models

# Extract the coordinates of selected atoms from all models
def extract_coordinates(models, atom_indices):
    coordinates = {index: [] for index in atom_indices}
    for model in models:
        for line in model:
            parts = line.split()
            atom_index = int(parts[1])  # Get atom index
            if atom_index in atom_indices:
                x = float(parts[5])  # X coordinate (6th column)
                y = float(parts[6])  # Y coordinate (7th column)
                z = float(parts[7])  # Z coordinate (8th column)
                coordinates[atom_index].append([x, y, z])
    return coordinates

# Choose which particle indices to export
def choose_particles(models):
    particle_indices = []
    print("Please choose the particle indices to export (separate with spaces):")
    print("You can find the particle indices in the corresponding .gro file.")
    
    # This part prints a list of particles available for selection, but the list won't be printed in the command line.
    input_indices = input("Enter indices: ").split()
    particle_indices = [int(idx) for idx in input_indices]
    return particle_indices

# Plot the 3D coordinates and save the plot as a PNG file
def plot_coordinates(coordinates, particle_index, output_dir):
    coords = np.array(coordinates[particle_index])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], label=f'Atom {particle_index}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title(f"3D Coordinates of Particle {particle_index}")
    
    # Save as PNG
    plt.savefig(f"{output_dir}/{particle_index}.png")
    plt.close()

# Plot the overlay of all particle coordinates and save the plot as a PNG file
def plot_all_coordinates(all_coordinates, output_dir):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for particle_index, coords in all_coordinates.items():
        coords = np.array(coords)
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], label=f'Atom {particle_index}')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title("Overlay of All Particles' Coordinates")
    
    # Save as PNG
    plt.savefig(f"{output_dir}/all_particles.png")
    plt.close()

# Save the coordinates of each particle to a TXT file
def save_coordinates_to_txt(coordinates, particle_index, output_dir):
    coords = np.array(coordinates[particle_index])
    np.savetxt(f"{output_dir}/{particle_index}.txt", coords, header="X Y Z", comments='')

# Main function to drive the program
def main():
    # Set input file paths and parameters
    gro_file = "npt-400.gro"  # Modify with your gro file path
    pdb_file = "5zncl-traj.pdb"  # Modify with your pdb file path
    output_dir = "output"  # Output folder

    # Read PDB file
    models = read_pdb(pdb_file)
    
    # Choose the particles to export
    particle_indices = choose_particles(models)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract coordinates of selected particles
    coordinates = extract_coordinates(models, particle_indices)
    
    # Save the coordinates and plot for each particle
    for particle_index in particle_indices:
        save_coordinates_to_txt(coordinates, particle_index, output_dir)
        plot_coordinates(coordinates, particle_index, output_dir)
    
    # Plot the overlay of all particles' coordinates
    plot_all_coordinates(coordinates, output_dir)
    
    print("Processing completed. Results have been saved to the output folder.")

if __name__ == "__main__":
    main()
