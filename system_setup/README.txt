# System setup scripts

This directory contains scripts and input files used to **build silica–water interface systems** for different simulations performed in the study *“Tailoring Heat Transfer at Silica–Water Interfaces via Hydroxyl and Methyl Surface Groups”* (V. Mandrolko et al., 2025).

---

## 📁 Contents

### 🔸 `SiO2_alpha_orthogonal_cell.data`
Initial atomic structure of crystalline α-SiO₂ used as the base unit for all generated systems.

### 🔸 `SiO2_base_for_ITR.py`
Python script that builds the silica slab (the base substrate) from the input data file.  
This structure serves as the foundation for constructing systems used in **interfacial thermal resistance (ITR)** simulations.

### 🔸 `water_slab_for_ITR.py`
Script that generates a slab of liquid water with user-defined thickness, used to couple with the silica base.

### 🔸 `system_for_ITR.py`
Combines the silica base and the water slab, applies the chosen surface functionalization (–OH / –CH₃), and produces the final system for ITR simulations.

### 🔸 Scripts for **adhesion** and **wetting angle**
Analogous scripts are provided for constructing systems used to compute the **work of adhesion** and the **contact angle (WA)**.  
The logic of the construction remains the same as for the ITR case, but with geometry and box dimensions adjusted for each specific simulation.

### 🔸 `forcefield_water_SiO2.potential`
Potential file used in the simulations. The same force field was applied across all constructed systems to ensure consistency.

---

## ⚙️ Usage notes

These scripts were written primarily for internal use and were not optimized for public release.  
They are provided *“as is”* to ensure transparency and reproducibility of the reported results.
