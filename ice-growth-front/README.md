# Ice Growth Front Detection

This folder contains two Python scripts used to analyze molecular dynamics simulation data to identify the ice growth front by analyzing particle distributions in layered slices along the z-axis.

These scripts design and methodology were inspired by:
Carignano, M. A., Shepson, P. B., & Szleifer, I. Molecular dynamics simulations of ice growth from supercooled water. *Molecular Physics*, **103**, 2957–2967 (2005).

## Overview

- `z_density_layer_analysis.py`
  Processes a single snapshot (GRO file) to calculate particle number distributions along the z-axis, dividing the simulation box into multiple layers. It outputs the layer-wise particle counts and identifies approximate boundaries between lower water, ice, and upper water regions based on the positions of `fSOL` residues.
- `ice_growth_front_detection.py`
  Processes a full trajectory (TRR + GRO files) frame-by-frame to detect the ice growth front dynamically. Using particle density deviations (with a given standard deviation from the first script), it identifies layers that belong to the ice phase and outputs the lower and upper ice layer boundaries over time.

## Recommended Setup

- A GRO file representing a snapshot configuration for density distribution analysis.
- A TRR trajectory file and matching GRO topology for the ice front detection across frames.

## Usage Instructions

### 1. `z_density_layer_analysis.py`

- Input: A single GRO file snapshot of your simulation system.
- Output:
  - `particle_distribution.txt`: Layer-wise particle counts.
  - `stat-result-total.txt`: Statistical summaries including mean and standard deviation.

Run the script directly:

```bash
python z_density_layer_analysis.py
```

### 2. `ice_growth_front_detection.py`

- Inputs:
  - TRR trajectory file (`npt-200.trr`)
  - GRO topology file (`npt-200.gro`)
- Output:
  - `1-sigma.txt`: Time series data containing the detected lower and upper ice layer boundaries with corresponding z-coordinates.

Before running, update parameters such as `num_layers`, `box_height`, and `std_dev` (obtained from the first script’s output).
Run the script:

```bash
python ice_growth_front_detection.py
```

## How the scripts work together

1. The first script analyzes a single snapshot to compute particle distributions per layer and calculate key statistics (mean, variance, standard deviation).
2. The computed standard deviation is used as a threshold in the second script.
3. The second script loads the full trajectory, calculates particle density per layer for each frame, and identifies ice layers based on deviations from the mean density using the threshold from the first script.
4. The output of the second script tracks the dynamic position of the ice front over simulation time.

## Test Files
Two example files, `npt-10.gro` and `npt-10.trr`, are provided for users to test and verify the basic functionality of the scripts. These files logically contain the necessary data to complete the intended analysis steps. However, due to GitHub's file size limitations, they have not been fully tested on real simulation data.

If you encounter any issues or unexpected behavior while using these test files or the scripts, please feel free to submit an issue on GitHub or contact via email for support.

## Note on Methodology
This script was independently developed based on literature insights and personal implementation during the early stages of the project. While it offers a novel approach to analyzing ice growth fronts, users are strongly encouraged to complement their analysis with widely validated mainstream tools such as CHILL+, which have undergone extensive benchmarking and community validation.
