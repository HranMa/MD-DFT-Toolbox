#!/usr/bin/env python
# coding: utf-8
'''
README
STEP 1: Create a new index.ndx file that includes all oxygen atoms from both WATER (SOL) and ICE (fSOL).
STEP 2: Generate a trajectory file that contains only the oxygen atoms. Use the new index file by specifying *-n index.ndx*. 
        Let's name the output file O-only.xtc. (traj.xtc and topol.tpr are examples, follow your setting)
        command: gmx trjconv -f traj.xtc -s topol.tpr -n index.ndx -o O-only.xtc
STEP 3: Analyze the O-only.xtc file by running the following command:
        python3 chill+_initial.py O-only.xtc

Details for Step 1:
1. Run the command: gmx make_ndx -f NPT.tpr -o index.ndx
2. Use the following commands in the interactive prompt to create two groups:
     - r SOL & a OW       → creates a group named "SOL_OW"
     - r fSOL & a fOW      → creates a group named "fSOL_OW"
   (These might be assigned to group numbers 10 and 11, for example.)
3. Merge the two groups by entering: 10 | 11
4. Rename the new group (e.g., group 12) to something meaningful like "fSOL_SOL_O".
5. Type "q" to save and quit.

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

