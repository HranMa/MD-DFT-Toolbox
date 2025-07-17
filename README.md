# MD & DFT Data Processing Scripts

This repository contains a collection of auxiliary scripts for data processing in Molecular Dynamics (MD) simulations and Density Functional Theory (DFT) calculations. These scripts were developed primarily during research under Prof. Wang Dan's guidance, with assistance from ChatGPT.

If you find these scripts useful in your research, please kindly cite the corresponding references under each section and star this repository to support the project.

---

## Project Structure

- **MSD**  
  Scripts for analyzing ion positions and mean squared displacement in simulations.  
  *Note:* It is generally recommended to use the built-in `gmx msd` tool from GROMACS for MSD calculations instead of these scripts.

- **chill+**  
  This script performs structural classification of water molecules in molecular dynamics trajectories using the CHILL+ algorithm, implemented via the OVITO Python API  

- **ice-growth-front**  
  Scripts to analyze GROMACS simulation data to determine ice growth fronts in systems where ice is centrally located, flanked by water layers.  
  *Details:*  
  - Residue names: ice = `fSOL`, water = `SOL`  
  - TIP4P/ice water model used; water molecules have 4 atoms in `.gro` files.

- **normal-tools**  
  A toolkit of small utility scripts for various routine tasks.

---

## Getting Started

### Environment Setup

- Recommended environment: WSL Ubuntu 22.04 or similar Linux distribution  
- Required dependencies: please check the import statements at the top of each Python script and install the corresponding packages accordingly.

---

## Contributing

Contributions and suggestions are welcome! 

Please open an issue or submit a pull request if you have improvements or bug fixes.

---

## Acknowledgements

Scripts were developed under Prof. Wang Dan’s supervision, with support from ChatGPT.

---

## Contact

For any questions or feedback, please contact: mahaorancn@163.com 

***Note**: I prefer communicating in Chinese, but English is totally fine.*
