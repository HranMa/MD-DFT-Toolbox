# normal-tools

This folder contains a set of lightweight Python scripts for structural or file-level preprocessing in molecular dynamics simulations.

## Scripts Overview

### 🧊 delete-free-water.py

Removes water molecules inside the ice seed zone or outside the simulation box in a `.gro` file.

- Recommended to match the box dimensions with the seed in at least two directions for periodicity.
- Do not overly rely on this script; prefer generating a clean initial structure.

### 🔁 hex-cubic-trans.py

Provides conversion between three-axis `(hkl)` and four-axis `(hkil)` crystallographic indices.

- Supports both plane and direction conversions.
- Includes interactive prompts for user input.

---

### 💤 Legacy / Reference Scripts

These are kept for archival purposes. Not actively used or maintained.

| Script                            | Description                                                  |
| --------------------------------- | ------------------------------------------------------------ |
| `zz_LEGACY_partial-stat.py`       | Early version for zone-wise particle statistics. Use newer tools if available. |
| `zz_LEGACY_random-coord.py`       | Generates random coordinates in `.gro`-compatible format. Reference only. |
| `zz_LEGACY_xvg-downsample-by2.py` | Halves the number of lines in a text-based trajectory to match output frequency. Useful only in specific mismatched sampling cases. |

---

## Notes

- Scripts are minimal and purpose-specific.
- No external dependencies beyond Python standard library and NumPy (for `partial-stat`).
- Always validate output manually if using any legacy tools.
