import random

"""
Note: This script may not be as effective as ChatGPT. Provided for reference only.

This script generates a specified number of random coordinates.
The output format follows the GROMACS .gro file style.
"""

def generate_random_coordinates(x_range, y_range, z_range, num_points=10, output_file="random-coord.txt"):
    # Open file for writing
    with open(output_file, "w") as file:
        for _ in range(num_points):
            # Generate three coordinate values within the given range
            x = round(random.uniform(x_range[0], x_range[1]), 3)
            y = round(random.uniform(y_range[0], y_range[1]), 3)
            z = round(random.uniform(z_range[0], z_range[1]), 3)
            # Write coordinates with 3 decimal places
            file.write(f"{x:.3f}   {y:.3f}   {z:.3f}\n")

# Hardcoded coordinate ranges (in nm)
x_range = (0, 3.54623)  # X range
y_range = (0, 3.54623)  # Y range
z_range = (3.35, 8.5)   # Z range

# Hardcoded number of output coordinate points
num_points = 15  # Number of random coordinates to generate

# Generate and write random coordinates
generate_random_coordinates(x_range, y_range, z_range, num_points, output_file="random-coord.txt")

print("Random coordinates have been generated and saved to 'random-coord.txt'.")
