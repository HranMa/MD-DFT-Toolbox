# CHILL+ Ice Structure Detection

This script performs structural classification of water molecules in molecular dynamics trajectories using the **CHILL+ algorithm**, implemented via the OVITO Python API. It identifies and tracks ice formation over time by selecting only oxygen atoms and analyzing the largest ice cluster in each frame.

> 📄 **References**:
>
> - **Nguyen, A. H. & Molinero, V.** *Identification of clathrate hydrates, hexagonal ice, cubic ice, and liquid water in simulations: The CHILL+ algorithm*. *J. Phys. Chem. B* **119**, 9369–9376 (2015).[https://doi.org/10.1021/acs.jpcb.5b04252](https://pubs.acs.org/doi/10.1021/jp510289t)
> - **Stukowski, A.** *Visualization and analysis of atomistic simulation data with OVITO–the Open Visualization Tool*. *Modelling Simul. Mater. Sci. Eng.* **18**, 015012 (2010). [https://doi.org/10.1088/0965-0393/18/1/015012](https://iopscience.iop.org/article/10.1088/0965-0393/18/1/015012)

------

## 🧊 Overview

This script:

- Uses the CHILL+ algorithm to classify each water molecule (oxygen atom) as ice or not.
- Filters out non-ice molecules.
- Applies cluster analysis to identify the largest contiguous ice structure.
- Outputs the size of the largest ice cluster frame-by-frame.

------

## 📦 Requirements

- Python ≥ 3.7

- OVITO Python module
   Install it via:

  ```bash
  pip install ovito
  ```

  For more information and API usage, refer to the [OVITO Python Reference — OVITO Python Reference 3.12.4 documentation](https://www.ovito.org/manual/python/)

------

## ▶️ Usage Instructions

### Step 1: Create an index file including oxygen atoms

The steps using `gmx make_ndx` are intended to **select all oxygen atoms** in the simulation box — including those in the **ice seed (`fSOL`)** and **bulk water (`SOL`)** — and merge them into a single index group. This unified group is required for CHILL+ to classify the local structure of each oxygen atom.

Use `gmx make_ndx` to group oxygen atoms from both water (`SOL`) and ice (`fSOL`):

```bash
gmx make_ndx -f topol.tpr -o index.ndx
```

In the interactive prompt:

```bash
r SOL & a OW
r fSOL & a fOW
10 | 11
name 12 fSOL_SOL_O
q
```

### Step 2: Extract oxygen-only trajectory

```bash
gmx trjconv -f traj.xtc -s topol.tpr -n index.ndx -o O-only.xtc
```

Replace filenames according to your project.

### Step 3: Run the script

```bash
python chill+.py O-only.xtc
```

------

## 📁 Output

- `ice_count.txt`: A tab-separated file with:

  | Frame | Timestep | Ice_Count |
  | ----- | -------- | --------- |

  

  Each line corresponds to the number of oxygen atoms in the **largest ice cluster** per frame.

------

## 🧪 Test File Provided

A sample file `O-only.xtc` is included in this folder to help users verify that the script is functioning correctly.

> ⚠️ **Note**:
>  While this sample data logically supports the workflow, it has not been fully validated due to GitHub’s file size limitations.
>  If you encounter issues or unexpected results, please open an issue or contact the author via email.

------

## 📖 Citations

If you use this code or workflow in your work, please cite both the CHILL+ algorithm and OVITO:

- Nguyen, A. H. & Molinero, V. *Identification of clathrate hydrates, hexagonal ice, cubic ice, and liquid water in simulations: The CHILL+ algorithm.*  
  *Journal of Physical Chemistry B*, **119**, 9369–9376 (2015).  
  [https://doi.org/10.1021/acs.jpcb.5b04252](https://pubs.acs.org/doi/10.1021/jp510289t)

- Stukowski, A. *Visualization and analysis of atomistic simulation data with OVITO—the Open Visualization Tool.*  
  *Modelling and Simulation in Materials Science and Engineering*, **18**, 015012 (2010).  
  [https://doi.org/10.1088/0965-0393/18/1/015012](https://iopscience.iop.org/article/10.1088/0965-0393/18/1/015012)
