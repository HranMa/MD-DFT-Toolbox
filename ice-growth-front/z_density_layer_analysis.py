import numpy as np

'''
This script is designed to work with ice_growth_front_detection.py.
It analyzes the ice growth front by computing the particle distribution along the z-axis in layers.
'''

def read_gro_file(gro_file):
    """Read a GRO file and return atom coordinates, box size,
    and min/max z-coordinates of 'fSOL' residues."""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    coordinates = []
    fSOL_z_coordinates = []
    
    # Read box dimensions from last line
    box_size = np.array([float(val) for val in lines[-1].split()])
    
    # Parse atom coordinates from the file
    for line in lines[2:]:
        if line.strip():
            try:
                res_name = line[5:10].strip()  # Residue name
                x = float(line[20:28].strip())
                y = float(line[28:36].strip())
                z = float(line[36:44].strip())
                coordinates.append([x, y, z])
                
                if res_name == 'fSOL':
                    fSOL_z_coordinates.append(z)
            except ValueError:
                print(f"Warning: Skipping line due to invalid coordinate data: {line.strip()}")
    
    z_fSOL_min = np.min(fSOL_z_coordinates) if fSOL_z_coordinates else None
    z_fSOL_max = np.max(fSOL_z_coordinates) if fSOL_z_coordinates else None
    
    return np.array(coordinates), box_size, z_fSOL_min, z_fSOL_max

def calculate_particle_distribution(coordinates, box_size, n_layers, z_fSOL_min, z_fSOL_max, output_file='particle_distribution.txt'):
    """Calculate particle counts per layer along z-axis and divide into water/ice regions."""
    z_min, z_max = np.min(coordinates[:, 2]), np.max(coordinates[:, 2])
    z_step = (z_max - z_min) / n_layers
    
    particle_distribution = np.zeros(n_layers)
    z_coordinates = np.zeros(n_layers)
    
    # Assign particles to layers by z coordinate
    for coord in coordinates:
        z = coord[2]
        layer_index = int((z - z_min) // z_step)
        if layer_index >= n_layers:
            layer_index = n_layers - 1
        particle_distribution[layer_index] += 1
    
    # Compute z coordinate representative of each layer
    for i in range(n_layers):
        if i == 0:
            z_coordinates[i] = z_min
        elif i == n_layers - 1:
            z_coordinates[i] = z_max
        else:
            z_coordinates[i] = z_min + (i + 0.5) * z_step
    
    # Define water and ice regions based on fSOL z-range
    lower_water_end = int((z_fSOL_min - z_min) // z_step)
    ice_layers_start = lower_water_end + 1
    ice_layers_end = int(ice_layers_start + (z_fSOL_max - z_fSOL_min) // z_step)
    upper_water_start = ice_layers_end + 1
    
    # Output layer particle distribution
    with open(output_file, 'w') as txtfile:
        txtfile.write("Layer Number\tz Coordinate (nm)\tParticle Count\n")
        for i in range(n_layers):
            txtfile.write(f"{i + 1}\t{z_coordinates[i]:.4f}\t{particle_distribution[i]}\n")
    
    print(f"\nData saved to '{output_file}'")
    print(f"\nLower Water Layers: 1 ~ {lower_water_end}")
    print(f"Ice Layers: {ice_layers_start} ~ {ice_layers_end}")
    print(f"Upper Water Layers: {upper_water_start} ~ {n_layers}")
    
    return lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers

def calculate_statistics(particle_distribution):
    """Calculate mean, variance, and standard deviation of particle counts."""
    mean = np.mean(particle_distribution)
    variance = np.var(particle_distribution)
    std_dev = np.std(particle_distribution)
    return mean, variance, std_dev

def output_statistics(mean, variance, std_dev, lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers, output_file="stat-result.txt"):
    """Save statistics and layer region info to a text file."""
    with open(output_file, 'w') as f:
        f.write(f"Mean: {mean:.4f}\n")
        f.write(f"Variance: {variance:.4f}\n")
        f.write(f"Standard Deviation: {std_dev:.4f}\n\n")
        f.write(f"Lower Water Layers: 1,{lower_water_end}\n")
        f.write(f"Ice Layers: {ice_layers_start},{ice_layers_end}\n")
        f.write(f"Upper Water Layers: {upper_water_start},{n_layers}\n")
    print(f"Statistics saved to '{output_file}'\n")

def main():
    """Main program entry point."""
    gro_file = 'npt-10.gro'  # Input GRO file path
    output_file = 'particle_distribution.txt'  # Output distribution file
    stat_file = 'stat-result-total.txt'  # Output statistics file
    n_layers = 200  # Number of layers
    
    coordinates, box_size, z_fSOL_min, z_fSOL_max = read_gro_file(gro_file)
    
    if z_fSOL_min is None or z_fSOL_max is None:
        print("Error: No 'fSOL' residues found in the GRO file.")
        return
    
    lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers = calculate_particle_distribution(
        coordinates, box_size, n_layers, z_fSOL_min, z_fSOL_max, output_file)
    
    particle_distribution = np.loadtxt(output_file, delimiter='\t', skiprows=1, usecols=2)
    mean, variance, std_dev = calculate_statistics(particle_distribution)
    
    print(f"\nMean Particle Count: {mean:.4f}")
    print(f"Standard Deviation: {std_dev:.4f}")
    
    output_statistics(mean, variance, std_dev, lower_water_end, ice_layers_start, ice_layers_end, upper_water_start, n_layers, stat_file)

if __name__ == "__main__":
    main()
