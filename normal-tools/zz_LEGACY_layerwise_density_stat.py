"""
DEPRECATED: This script is outdated and rarely used. Kept only for reference or backup.

Note:
This script does NOT overwrite the output file by default.
If executed multiple times, duplicate statistics will accumulate. 
Please manually remove unwanted content from the output file if needed.

The current logic appends statistics in the following order:
Upper Water → Ice → Lower Water
"""

import numpy as np

def read_particle_distribution(file_path):
    """Read particle count from input file (npt-400.txt); return per-layer counts."""
    layers = []
    particle_counts = []

    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip header

        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                layers.append(int(parts[0]))
                try:
                    particle_count = float(parts[2])
                    if particle_count >= 0:
                        particle_counts.append(particle_count)
                    else:
                        print(f"Warning: Negative particle count at layer {parts[0]}, skipping.")
                except ValueError:
                    print(f"Warning: Invalid particle count at layer {parts[0]}, skipping.")

    return np.array(layers), np.array(particle_counts)

def read_layer_info(stat_file):
    """Read layer range definitions from stat-result-total.txt."""
    layer_info = {}

    with open(stat_file, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if "Lower Water Layers" in line:
                lower_range = line.split(":")[-1].strip()
                layer_info['lower_water'] = list(map(int, lower_range.split(',')))
            elif "Ice Layers" in line:
                ice_range = line.split(":")[-1].strip()
                layer_info['ice'] = list(map(int, ice_range.split(',')))
            elif "Upper Water Layers" in line:
                upper_range = line.split(":")[-1].strip()
                layer_info['upper_water'] = list(map(int, upper_range.split(',')))

    return layer_info

def calculate_statistics(particle_counts):
    """Compute mean, variance, and standard deviation."""
    mean = np.mean(particle_counts)
    variance = np.var(particle_counts)
    std_dev = np.std(particle_counts)
    return mean, variance, std_dev

def output_statistics(mean, variance, std_dev, region_name, output_file):
    """Append computed statistics for a region to output file."""
    with open(output_file, 'a') as f:
        f.write(f"{region_name} Mean: {mean:.4f}\n")
        f.write(f"{region_name} Variance: {variance:.4f}\n")
        f.write(f"{region_name} Standard Deviation: {std_dev:.4f}\n\n")

def main():
    # Input files
    npt_file = 'npt-400.txt'               # Output from z-density-stat-partial.py
    stat_file = 'stat-result-total.txt'    # Layer range info from z-density-stat-partial.py

    # Output file
    output_file = 'partial-stat.txt'

    layers, particle_counts = read_particle_distribution(npt_file)
    layer_info = read_layer_info(stat_file)

    lower_range = layer_info.get('lower_water', [])
    ice_range = layer_info.get('ice', [])
    upper_range = layer_info.get('upper_water', [])

    # Clamp ranges to array length
    if lower_range and lower_range[1] > len(particle_counts):
        lower_range[1] = len(particle_counts)
    if ice_range and ice_range[1] > len(particle_counts):
        ice_range[1] = len(particle_counts)
    if upper_range and upper_range[1] > len(particle_counts):
        upper_range[1] = len(particle_counts)

    # Extract per-region counts
    lower_counts = particle_counts[lower_range[0]-1 : lower_range[1]] if lower_range else []
    ice_counts = particle_counts[ice_range[0]-1 : ice_range[1]] if ice_range else []
    upper_counts = particle_counts[upper_range[0]-1 : upper_range[1]] if upper_range else []

    # Output statistics
    if lower_counts.size > 0:
        mean, var, std = calculate_statistics(lower_counts)
        output_statistics(mean, var, std, 'Lower Water', output_file)

    if ice_counts.size > 0:
        mean, var, std = calculate_statistics(ice_counts)
        output_statistics(mean, var, std, 'Ice', output_file)

    if upper_counts.size > 0:
        mean, var, std = calculate_statistics(upper_counts)
        output_statistics(mean, var, std, 'Upper Water', output_file)

    print(f"Partial statistics saved to '{output_file}'")

if __name__ == "__main__":
    main()
