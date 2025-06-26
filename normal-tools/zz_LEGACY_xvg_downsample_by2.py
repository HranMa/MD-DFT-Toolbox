"""
IMPORTANT:
- Do NOT overly rely on this script for resampling time series data.
- Best practice is to ensure consistent .mdp settings (e.g., nstxout, nstenergy) before running MD simulations.

Script: xvg_downsample_by2.py
Purpose:
    This script removes every other line from a text-based time series file 
    (e.g., GROMACS .xvg or exported .txt), effectively halving the sampling frequency.

Usage:
    - Use this when one trajectory was saved at double the frequency (e.g., 0.5 ns/frame),
      while others were saved at 1.0 ns/frame.
    - It keeps every second line (i.e., odd-numbered lines), starting from line 0.

Input:
    input_file  = name of the input file (e.g., 'water-ice.txt')
    output_file = name of the output file (e.g., 'water-ice-output.txt')

Note:
    Only works for text files with one data point per line, and assumes the first line is to be kept.
"""

input_file = 'water-ice.txt'  # Replace with your actual input file name
output_file = 'water-ice-output.txt'  # Output file name

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for index, line in enumerate(infile):
        # Write only even-indexed lines (0-based), i.e., every other line
        if index % 2 == 0:
            outfile.write(line)

print(f"Even-indexed lines retained. Output saved to '{output_file}'.")
