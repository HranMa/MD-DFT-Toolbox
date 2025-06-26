import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI Agg backend for matplotlib
import matplotlib.pyplot as plt

'''
NOTE:
- This script is designed to be used together with particle-position.py.
- However, it is generally NOT recommended to use these scripts for MSD calculation.
- Instead, use the built-in Gromacs tool 'gmx msd' for accurate and efficient analysis.
'''

# Set box dimensions for periodic boundary conditions (Lx, Ly, Lz) in Angstroms
Lx, Ly, Lz = 35, 35, 90

# Function to apply periodic boundary correction
def apply_periodic_boundary(x, L):
    return x - L * np.round(x / L)

# Read particle trajectory data from a file
def read_particle_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)  # Skip header line
    return data  # Returns an (n,3) array where n=time steps, 3=X,Y,Z coords

# Calculate Mean Squared Displacement (MSD)
def calculate_msd(trajectory):
    initial_position = trajectory[0]  # Position at time zero
    displacements = trajectory - initial_position  # Displacements at all times
    # Apply periodic boundary correction
    displacements[:, 0] = apply_periodic_boundary(displacements[:, 0], Lx)
    displacements[:, 1] = apply_periodic_boundary(displacements[:, 1], Ly)
    displacements[:, 2] = apply_periodic_boundary(displacements[:, 2], Lz)

    msd = np.mean(np.square(displacements), axis=1)
    return msd

# Get all txt files from output directory
def get_txt_files():
    files = glob.glob("./output/*.txt")
    return files

def main():
    txt_files = get_txt_files()

    print("The following files will be processed:")
    for f in txt_files:
        print(f)
    
    proceed = input("Do you want to continue processing these files? (y/n): ")
    if proceed.lower() != 'y':
        print("Processing canceled.")
        return

    msd_all = []  # Store MSD arrays for all files

    for file_path in txt_files:
        print(f"Processing file: {file_path}")
        trajectory = read_particle_data(file_path)
        msd = calculate_msd(trajectory)
        msd_all.append(msd)

    msd_all = np.array(msd_all)
    mean_msd = np.mean(msd_all, axis=0)

    os.makedirs('./msd', exist_ok=True)
    np.savetxt('./msd/msd.txt', mean_msd)

    plt.figure()
    plt.plot(mean_msd, label="MSD")
    plt.xlabel("Time Steps")
    plt.ylabel("Mean Squared Displacement (MSD) [Å²]")
    plt.title("Mean Squared Displacement (MSD)")
    plt.legend()
    plt.savefig('./msd/msd.png', dpi=300)

    print("Processing complete. MSD results saved to msd/msd.txt and msd/msd.png")

if __name__ == "__main__":
    main()
