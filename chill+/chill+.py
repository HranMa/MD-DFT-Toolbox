#!/usr/bin/env python
# coding: utf-8
'''
CHILL+ Ice Structure Analysis Script

This script identifies and quantifies ice structures from a molecular dynamics trajectory
using the Chill+ algorithm implemented in OVITO.

Usage Overview:

STEP 1: Create an index file including all oxygen atoms from both WATER (SOL) and ICE (fSOL).
    Example:
    1. Run:
        gmx make_ndx -f NPT.tpr -o index.ndx
    2. In the interactive prompt, create groups for oxygen atoms:
        - r SOL & a OW      (group name e.g. "SOL_OW")
        - r fSOL & a fOW    (group name e.g. "fSOL_OW")
    3. Merge these two groups (e.g., "10 | 11")
    4. Rename the merged group (e.g., "fSOL_SOL_O")
    5. Save and quit (type "q").

STEP 2: Extract a trajectory containing only oxygen atoms with the new index file.
    Example command:
        gmx trjconv -f traj.xtc -s topol.tpr -n index.ndx -o O-only.xtc

STEP 3: Run this script on the oxygen-only trajectory:
        python3 chill+_initial.py O-only.xtc

Output:
- The script writes 'ice_count.txt', containing frame-by-frame ice particle counts of the largest ice cluster.
- Console output provides progress with frame number, timestep, and ice count.

Dependencies:
- OVITO Python API
- GROMACS tools (for indexing and trajectory extraction)

'''

from ovito.io import import_file
from ovito.modifiers import (
    ChillPlusModifier, SelectTypeModifier, DeleteSelectedModifier,
    ClusterAnalysisModifier, ExpressionSelectionModifier
)
import sys

# Get input trajectory file name from command-line argument
trjfilename = sys.argv[1]

# Output file name (fixed to txt)
output_filename = 'ice_count.txt'

# Load the trajectory
pipeline = import_file(trjfilename)

# Add Chill+ modifier to identify ice structures
pipeline.modifiers.append(ChillPlusModifier())

# Remove non-ice particles (0: Other, 4: Unknown, 5: Liquid)
pipeline.modifiers.append(SelectTypeModifier(
    property="Structure Type",
    types={0, 4, 5}
))
pipeline.modifiers.append(DeleteSelectedModifier())

# Perform cluster analysis to detect ice clusters
pipeline.modifiers.append(ClusterAnalysisModifier(
    cutoff=2.9,
    sort_by_size=True,
    compute_com=True
))

# Keep only the largest ice cluster (Cluster ID = 1)
pipeline.modifiers.append(ExpressionSelectionModifier(expression="Cluster != 1"))
pipeline.modifiers.append(DeleteSelectedModifier())

# Write results to a txt file (Tab-separated, Origin-friendly)
with open(output_filename, 'w') as file_out:
    file_out.write('Frame\tTimestep\tIce_Count\n')

    for nframe in range(pipeline.source.num_frames):
        data = pipeline.compute(nframe)
        timestep = data.attributes['Timestep']
        num_particles = data.particles.count
        print(f'Frame {nframe}: Timestep = {timestep}, Ice_Count = {num_particles}')
        file_out.write(f'{nframe}\t{timestep:.1f}\t{num_particles}\n')
