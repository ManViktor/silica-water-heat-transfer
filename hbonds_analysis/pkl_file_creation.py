# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 20:59:46 2025

@author: mandrolk1
"""

import os, math, pickle
import numpy as np
from scipy.spatial import cKDTree

# Parameters
dist_cutoff = 2.5
angle_cutoff = 150.0
n_steps = 10000
box = np.array([60.333384, 139.373984, 1000.0])
base_dir = os.path.dirname(os.path.abspath(__file__))  # Use the folder where the script is located
percentages = ['70_perc']

def load_pairs(filename):
    if not os.path.exists(filename):
        return np.empty((0, 2), dtype=int)
    return np.loadtxt(filename, dtype=int, skiprows=1)

def load_water_molecules(filename):
    if not os.path.exists(filename):
        return np.empty((0, 3), dtype=int)
    return np.loadtxt(filename, dtype=int, skiprows=1)

def parse_timesteps(xyz_file, n_steps):
    with open(xyz_file, 'r') as f:
        for _ in range(n_steps):
            if not f.readline().startswith('ITEM: TIMESTEP'):
                raise ValueError("Invalid xyz format")
            ts = int(f.readline().strip())
            f.readline(); natoms = int(f.readline().strip())
            f.readline()
            for _ in range(3): f.readline()
            f.readline()
            ids = np.empty(natoms, int)
            names = np.empty(natoms, '<U2')
            coords = np.empty((natoms, 3), float)
            for i in range(natoms):
                parts = f.readline().split()
                ids[i] = int(parts[0]); names[i] = parts[2]
                coords[i] = list(map(float, parts[3:6]))
            yield {'ts': ts, 'ids': ids, 'names': names, 'coords': coords}

def compute_angles_vec(coords1, coords2, coords3, box):
    v1 = coords1 - coords2
    v2 = coords3 - coords2
    for i in (0, 1):
        v1[:, i] -= box[i] * np.round(v1[:, i] / box[i])
        v2[:, i] -= box[i] * np.round(v2[:, i] / box[i])
    cosang = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    return np.degrees(np.arccos(np.clip(cosang, -1, 1)))

def get_slabs(frame, perc):
    ids, names, coords = frame['ids'], frame['names'], frame['coords']
    mask = (names == 'O2') if perc == '100_perc' else (names == 'C')
    z = coords[mask, 2]
    low = z[z < 80].mean() - 1.5; slab_low = (low, low + 5)
    high = z[z > 80].mean() - 1.5; slab_high = (high - 5, high)
    return slab_low, slab_high

def process_system(perc):
    print(f"Processing {perc}")
    res = os.path.join(base_dir)
    xyz = os.path.join(res, "dynamics.xyz")
    oh = load_pairs(os.path.join(res, "oh_groups.dat"))
    ch3 = load_pairs(os.path.join(res, "ch3_groups.dat"))
    water_mols = load_water_molecules(os.path.join(res, "water_molecules_ids.dat"))
    slab_low, slab_high = get_slabs(next(parse_timesteps(xyz, 1)), perc)

    hbonds_per_timestep = {}
    for frame in parse_timesteps(xyz, n_steps):
        ts, ids, names, coords = frame.values()
        z = coords[:, 2]
        mask_w = (names == 'O') & (((z >= slab_low[0]) & (z <= slab_low[1])) | ((z >= slab_high[0]) & (z <= slab_high[1])))
        ow_ids = ids[mask_w]; ow_coords = coords[mask_w]
        ow_xy = np.mod(ow_coords[:, :2], box[:2])
        tree_w = cKDTree(ow_xy, boxsize=box[:2].tolist())

        hbonds_OH_donor = set()
        hbonds_CH3_donor = set()
        hbonds_OH_acceptor = set()

        # OH donor
        for o2, h1 in oh:
            hi = np.where(ids == h1)[0]; oi = np.where(ids == o2)[0]
            if hi.size and oi.size:
                h_coord = coords[hi[0], :2] % box[:2]
                o2_coord = coords[oi[0]]
                for j in tree_w.query_ball_point(h_coord, dist_cutoff):
                    owp = ow_coords[j]
                    dx = h_coord - owp[:2]; dx -= box[:2] * np.round(dx / box[:2])
                    if np.hypot(*dx) <= dist_cutoff:
                        angles = compute_angles_vec(o2_coord[np.newaxis], coords[hi[0]][np.newaxis], owp[np.newaxis], box)
                        if angles[0] > angle_cutoff:
                            hbonds_OH_donor.add((int(h1), int(ids[mask_w][j])))

        # CH3 donor
        bondedC = set()
        for c, h2 in ch3:
            if c in bondedC: continue
            ci = np.where(ids == c)[0]; hi = np.where(ids == h2)[0]
            if ci.size and hi.size:
                c_coord = coords[ci[0]]; h2_coord = coords[hi[0], :2] % box[:2]
                for j in tree_w.query_ball_point(h2_coord, dist_cutoff):
                    owp = ow_coords[j]; dx = h2_coord - owp[:2]
                    dx -= box[:2] * np.round(dx / box[:2])
                    if np.hypot(*dx) <= dist_cutoff:
                        angles = compute_angles_vec(c_coord[np.newaxis], coords[hi[0]][np.newaxis], owp[np.newaxis], box)
                        if angles[0] > angle_cutoff:
                            hbonds_CH3_donor.add((int(h2), int(ids[mask_w][j])))
                            bondedC.add(c); break

        # OH acceptor
        water_sl = water_mols[np.isin(water_mols[:, 0], ow_ids)]
        h_ids = []; h_to_ow = {}
        for ow_id, h1_id, h2_id in water_sl:
            h_ids.extend([h1_id, h2_id])
            h_to_ow[h1_id] = ow_id; h_to_ow[h2_id] = ow_id

        if h_ids:
            h_coords = []
            for hid in h_ids:
                idx = np.where(ids == hid)[0]
                if idx.size:
                    h_coords.append((hid, coords[idx[0]]))
            h_xy = np.array([hc[1][:2] for hc in h_coords]) % box[:2]
            tree_h = cKDTree(h_xy, boxsize=box[:2].tolist())
            for o2, _ in oh:
                o2_idx = np.where(ids == o2)[0]
                if not o2_idx.size: continue
                o2_coord = coords[o2_idx[0]]
                o2_xy = o2_coord[:2] % box[:2]
                for j in tree_h.query_ball_point(o2_xy, dist_cutoff):
                    hid, h_full = h_coords[j]
                    ow_id = h_to_ow[hid]
                    ow_idx = np.where(ids == ow_id)[0]
                    if not ow_idx.size: continue
                    ow_full = coords[ow_idx[0]]
                    dx = o2_xy - (h_full[:2] % box[:2]); dx -= box[:2] * np.round(dx / box[:2])
                    if np.hypot(*dx) > dist_cutoff: continue
                    angles = compute_angles_vec(o2_coord[np.newaxis], h_full[np.newaxis], ow_full[np.newaxis], box)
                    if angles[0] > angle_cutoff:
                        hbonds_OH_acceptor.add((int(hid), int(o2)))

        hbonds_per_timestep[ts] = {
            "OH_donor": hbonds_OH_donor,
            "CH3_donor": hbonds_CH3_donor,
            "OH_acceptor": hbonds_OH_acceptor
        }
        print(f"Step {ts}: OH_donor={len(hbonds_OH_donor)}, CH3_donor={len(hbonds_CH3_donor)}, OH_acceptor={len(hbonds_OH_acceptor)}")

    with open(os.path.join(res, "hbonds_per_timestep.pkl"), "wb") as f:
        pickle.dump(hbonds_per_timestep, f)

if __name__ == '__main__':
    for perc in percentages:
        process_system(perc)
    print('All done.')