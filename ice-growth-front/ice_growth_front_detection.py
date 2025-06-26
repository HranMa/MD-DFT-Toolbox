import numpy as np
import mdtraj as md
from tqdm import tqdm  # Progress bar library

'''
This script works together with z_density_layer_analysis.py.
It processes the trajectory frame-by-frame to identify the ice growth front
by detecting the lower and upper ice layers based on particle density deviations.
'''

# Constants and input/output file paths
trr_file = 'npt-200.trr'  # Trajectory file path
top_file = 'npt-200.gro'  # GRO topology file path
output_file = '1-sigma.txt'  # Output results file path
num_layers = 200  # Number of layers to divide the simulation box into
box_height = 8.90420  # Total box height in nm (from GRO file)
std_dev = 34.0878  # Standard deviation from z_density_layer_analysis.py results

# Load trajectory using MDTraj
traj = md.load(trr_file, top=top_file)
time_steps = traj.time  # Simulation time points
z_coords = traj.xyz[:, :, 2]  # Extract z-coordinates of all atoms per frame
num_atoms = traj.n_atoms

def calculate_layer_density(z_coords, num_layers, box_height):
    """
    Calculate the particle count in each z-layer for a single frame.
    Args:
        z_coords: Array of z-coordinates for all atoms in the frame
        num_layers: Total number of layers dividing the box height
        box_height: Height of the simulation box in nm
    Returns:
        layer_density: Particle count per layer
        layer_z: Boundaries of each z-layer
    """
    layer_density = np.zeros(num_layers)
    layer_z = np.linspace(0, box_height, num_layers + 1)  # Layer boundaries along z
    
    for atom_z in z_coords:
        layer_index = int(atom_z / (box_height / num_layers))
        layer_index = min(layer_index, num_layers - 1)  # Ensure index is valid
        layer_density[layer_index] += 1

    return layer_density, layer_z

def define_ice_layers(layer_density, num_layers):
    """
    Identify which layers correspond to ice by detecting significant
    deviations in particle density compared to average.
    Args:
        layer_density: Particle counts per layer
        num_layers: Total number of layers
    Returns:
        ice_layers: Boolean array indicating ice layers
    """
    avg_density = np.mean(layer_density)
    ice_layers = np.zeros(num_layers, dtype=bool)
    
    for i in range(num_layers):
        if (layer_density[i] > (avg_density + np.sqrt(2) * std_dev)) or \
           (layer_density[i] < (avg_density - np.sqrt(2) * std_dev)):
            ice_layers[i] = True

    return ice_layers

def find_lower_ice_layer(ice_layers, midpoint):
    """
    Find the first ice layer from the bottom up to the midpoint.
    Args:
        ice_layers: Boolean array of ice layers
        midpoint: Middle layer index (geometric center)
    Returns:
        Index of the lower ice layer or -1 if none found
    """
    for i in range(midpoint):
        if ice_layers[i]:
            return i
    return -1

def find_upper_ice_layer(ice_layers, midpoint, num_layers, layer_density):
    """
    Find the first ice layer from near the top down to the midpoint.
    Args:
        ice_layers: Boolean array of ice layers
        midpoint: Middle layer index
        num_layers: Total number of layers
        layer_density: Particle counts per layer
    Returns:
        Index of the upper ice layer or midpoint if none found
    """
    avg_density = np.mean(layer_density)
    start_layer = max(num_layers - 20, 0)
    for i in range(start_layer, midpoint - 1, -1):
        if (layer_density[i] > (avg_density + np.sqrt(2) * std_dev)) or \
           (layer_density[i] < (avg_density - np.sqrt(2) * std_dev)):
            return i
    return midpoint

def process_trajectory(traj, num_layers, box_height, output_file):
    """
    Process all frames in the trajectory to identify ice layer boundaries
    and save the results to a file.
    Args:
        traj: MDTraj trajectory object
        num_layers: Number of layers dividing the box height
        box_height: Height of the simulation box in nm
        output_file: Path to save the results
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("Time (ns)\tLower Ice Layer\tUpper Ice Layer\tLower Z Coordinate (nm)\tUpper Z Coordinate (nm)\n")

        for step in tqdm(range(len(traj)), desc="Processing trajectory", unit="frame"):
            z_coords_step = z_coords[step]
            layer_density, layer_z = calculate_layer_density(z_coords_step, num_layers, box_height)
            ice_layers = define_ice_layers(layer_density, num_layers)
            midpoint = num_layers // 2
            
            lower_ice_layer = find_lower_ice_layer(ice_layers, midpoint)
            upper_ice_layer = find_upper_ice_layer(ice_layers, midpoint, num_layers, layer_density)

            # If no ice layers found, set defaults
            if lower_ice_layer == -1:
                lower_ice_layer = 0
            if upper_ice_layer == -1:
                upper_ice_layer = num_layers - 1

            lower_z = layer_z[lower_ice_layer] if lower_ice_layer != 0 else 0
            upper_z = layer_z[upper_ice_layer]
            time_ns = traj.time[step] / 1000.0  # Convert ps to ns

            f.write(f"{time_ns:.3f}\t{lower_ice_layer + 1}\t{upper_ice_layer + 1}\t{lower_z:.3f}\t{upper_z:.3f}\n")

# Run the processing function
process_trajectory(traj, num_layers, box_height, output_file)

print(f"Results have been saved to {output_file}")
