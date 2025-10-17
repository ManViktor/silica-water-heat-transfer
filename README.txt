This repository contains scripts and input files used in the study  
**"Tailoring Heat Transfer at Silica–Water Interfaces via Hydroxyl and Methyl Surface Groups"**  
(Viktor Mandrolko, Konstantinos Termentzidis, David Lacroix and Mykola Isaiev, 2025)

The repository is organized into two main parts.

---

## 📂 Repository structure

### `system_setup/`
This folder contains scripts related to the **construction of silica–water interface systems**.  
The scripts were used for generating the configurations employed in the evaluation of:

- **Contact angle** (wettability)
- **Work of adhesion**
- **Interfacial thermal resistance (ITR)**

In this folder, you will also find:
- the **potential file** used in the simulations,  
- and *silica structure file** used in the study.

> Note: these scripts were not originally written for public release, and may not be fully optimized or user-friendly. They are shared here “as is” for the sake of transparency and reproducibility.

---

### `hbonds_analysis/`
This folder contains analysis tools for **hydrogen-bond dynamics and interfacial structure**.  
The scripts were used to quantify hydrogen-bond lifetimes and formation/rupture rates at the functionalized silica surfaces.

---

## General information
- **Language:** Python 3  
- **Dependencies:** `numpy`, `matplotlib`, `pickle`, `collections`, `scipy`  

---

## 📚 Citation
If you use or adapt these scripts, please cite:

> Mandrolko, V., Termentzidis, K., Lacroix, D. and Isaiev, M. *Tailoring Heat Transfer at Silica–Water Interfaces via Hydroxyl and Methyl Surface Groups* (2025).  
> GitHub repository: [https://github.com/viktor-mandroliuk/silica-water-heat-transfer]

---
