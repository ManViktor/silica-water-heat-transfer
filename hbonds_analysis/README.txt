📁 Structure overview
hbonds_analysis/
│
├── find_water_molecules_ids.py
├── find_functional_molecules_ids.py
├── pkl_file_creation.py
├── continuous_ACF_lifetime.py
└── dynamics.xyz  ← trajectory file from LAMMPS

Step 1. Identify molecular groups

Scripts:

find_water_molecules_ids.py
find_functional_molecules_ids.py

These scripts parse the LAMMPS trajectory file (dynamics.xyz) to identify:

water molecules (O–H–H groups),
surface functional groups (–OH and –CH₃).

They generate text files containing molecule/atom IDs:

water_molecules_ids.dat
oh_groups.dat
ch3_groups.dat

These files are used later to locate hydrogen bonds between water and surface sites.

Step 2. Compute hydrogen bonds per timestep

Script: pkl_file_creation.py

This script:
Reads the .xyz trajectory and the ID files from Step 1.
Detects hydrogen bonds frame-by-frame based on:
Distance cutoff: 2.5 Å
Angle cutoff: 150°
Considers both donor and acceptor roles of functional groups and water.

Separately tracks OH_donor, OH_acceptor and CH3_donor

The output is a Pickle file storing hydrogen bond information for all timesteps:

hbonds_per_timestep.pkl

This file is the main input for the next analysis step.

Step 3. Compute hydrogen bond lifetime (Continuous ACF method)

Scripts: 

continuous_ACF_lifetime.py
intermittent_ACF_lifetime.py

This script performs a continuous/intermittent autocorrelation function (ACF) analysis on the hydrogen bond dynamics.

Main features:

Splits the trajectory into blocks and averages ACFs across them.
Calculates mean lifetime τ for each hydrogen bond type:
Fits an exponential decay to obtain τ, or uses integral fallback.

Produces:

.dat file with average ACFs
.png plots of ACF curves and exponential fits
console output with τ_mean ± τ_std

Adjustable parameters inside the script:

dt = 0.01          # ps between saved frames
block_size = 200   # number of steps per block
nsteps = 10000     # number of steps to process
num_processes = 1  # parallel cores to use
