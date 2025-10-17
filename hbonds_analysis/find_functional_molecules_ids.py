# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 11:23:20 2025

@author: mandrolk1
"""

# -*- coding: utf-8 -*-
"""
Adapted on Mon Jul  7 2025
Group functional hydroxyl (O2–H1) and methyl (C–H2) pairs from first timestep
Print atom counts and write indices of each pair to output files.
Improved PBC handling by replicating partner atoms in periodic images.
"""
import numpy as np
from scipy.spatial import cKDTree

def parse_first_timestep(input_filename):
    atoms = []
    with open(input_filename, 'r') as f:
        inside = False
        atoms_section = False
        for line in f:
            line = line.strip()
            if line.startswith("ITEM: TIMESTEP"):
                if not inside:
                    inside = True
                else:
                    break
                continue
            if inside:
                if line.startswith("ITEM: ATOMS"):
                    atoms_section = True
                    continue
                if atoms_section and line:
                    parts = line.split()
                    if len(parts) >= 6:
                        atoms.append((int(parts[0]), parts[2],
                                      float(parts[3]), float(parts[4]), float(parts[5])))
    dtype = [('id', 'i4'), ('name', 'U4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')]
    return np.array(atoms, dtype=dtype)

def compute_distance(a, b, box):
    delta = a - b
    for i in (0,1):
        delta[i] -= box[i] * np.round(delta[i] / box[i])
    return np.linalg.norm(delta)

def build_env_tree(partners, box):
    """Replicate partner coords in -1,0,+1 shifts for x,y; return KDTree and list of (orig_index)"""
    shifts = np.array([(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1)])
    coords = np.vstack((partners['x'], partners['y'], partners['z'])).T
    env_coords = []
    env_map = []
    for idx, coord in enumerate(coords):
        for shift in shifts:
            xy = coord[:2] + shift * box[:2]
            env_coords.append((xy[0], xy[1], coord[2]))
            env_map.append(idx)
    env_coords = np.array(env_coords)
    tree = cKDTree(env_coords[:,:2])
    return tree, env_coords, env_map

def find_pairs(atom_data, center_name, partner_name, cutoff, box, k_neighbors, output_file):
    centers = atom_data[atom_data['name']==center_name]
    partners = atom_data[atom_data['name']==partner_name]
    if center_name=='O2' and partner_name=='H1':
        total=len(atom_data)
        cnts={nm:np.sum(atom_data['name']==nm) for nm in('C','H2','O2','H1')}
        print(f"Total atoms: {total}")
        for nm,c in cnts.items(): print(f"{nm}: {c}")
    if len(centers)==0 or len(partners)==0:
        print(f" {center_name} {partner_name}")
        return
    tree, env_coords, env_map = build_env_tree(partners, box)
    results=[]
    for c in centers:
        c_xyz = np.array((c['x'],c['y'],c['z']))
        # query candidates within cutoff in xy
        idxs = tree.query_ball_point((c['x'],c['y']), r=cutoff)
        # compute full distances to original partners
        dlist=[]
        for env_i in idxs:
            p_i = env_map[env_i]
            p = partners[p_i]
            p_xyz = np.array((p['x'],p['y'],p['z']))
            r = compute_distance(c_xyz, p_xyz, box)
            if r<=cutoff:
                dlist.append((r, p_i))
        # sort by distance and take up to k_neighbors
        dlist.sort(key=lambda x: x[0])
        for _, p_i in dlist[:k_neighbors]:
            results.append((c['id'], partners[p_i]['id']))
    with open(output_file,'w') as out:
        out.write(f"{center_name}\t{partner_name}\n")
        for a,b in results: out.write(f"{a}\t{b}\n")
    print(f"{len(results)}  {center_name}-{partner_name},  {output_file}")

if __name__=='__main__':
    seed = 6
    perc = 50
    input_file = "dynamics.xyz"
    box = np.array([60.333384,139.333984,np.inf])
    cutoff=1.2
    find_pairs(parse_first_timestep(input_file),'O2','H1',cutoff,box,1,
               "oh_groups.dat")
    find_pairs(parse_first_timestep(input_file),'C','H2',cutoff,box,3,
               "ch3_groups.dat")
