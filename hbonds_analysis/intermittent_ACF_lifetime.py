# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 10:00:45 2025

@author: mandrolk1
"""

"""
Intermittent ACF (only) for hydrogen bonds.
- reads hbonds_per_timestep.pkl produced by your preprocessing
- builds per-pair presence arrays (0/1)
- computes intermittent ACF robustly and clamps to [0,1]
- computes integrated lifetime (tau = ∫ C(t) dt)
- optional: computes per-pair run-length lifetimes and histogram
- saves: .dat file with all C(t) columns, png plot, summary txt, lifetimes csv

Author: adapted for Viktor (Pomahайко)
"""
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


# ---------------- USER SETTINGS ----------------
dt = 0.01           # ps between saved frames
skip_factor = 1     
nsteps = 10000       
block_size = 1000   # number of steps per block
bond_types = ["OH_donor", "OH_acceptor", "CH3_donor"]
verbose = True
save_pair_lifetimes = True
# ------------------------------------------------
def build_presence_arrays(hbonds_per_timestep, timesteps, bt):
    """
    Return (arr, pairs_list)
    arr shape = (n_pairs, n_timesteps), dtype = int8 (0/1)
    pairs_list = sorted list of pairs (tuples)
    """
    pairs_set = set()
    for ts in timesteps:
        entry = hbonds_per_timestep[ts]
        if isinstance(entry, dict):
            current = entry.get(bt, set())
        elif isinstance(entry, set):
            current = entry if bt == "all" else set()
        else:
            raise TypeError(f"Unexpected type at timestep {ts}: {type(entry)}")
        pairs_set.update(current)

    pairs_list = sorted(pairs_set)
    n_pairs = len(pairs_list)
    n_ts = len(timesteps)
    if n_pairs == 0:
        return np.zeros((0, n_ts), dtype=np.int8), pairs_list

    arr = np.zeros((n_pairs, n_ts), dtype=np.int8)
    # fill rows
    for j, ts in enumerate(timesteps):
        entry = hbonds_per_timestep[ts]
        if isinstance(entry, dict):
            current = entry.get(bt, set())
        elif isinstance(entry, set):
            current = entry if bt == "all" else set()
        # convert current to a set for membership
        for i, pair in enumerate(pairs_list):
            if pair in current:
                arr[i, j] = 1
    return arr, pairs_list


def compute_intermittent_acf_from_arr(arr, verbose=False):
    """
    arr: shape (n_pairs, n_ts), dtype 0/1 (int8 or bool)
    Computes intermittent ACF:
      For each tau: C(tau) = sum_p sum_{t=0}^{N-1-tau} arr[p,t]*arr[p,t+tau] /
                         sum_p sum_{t=0}^{N-1-tau} arr[p,t]
    Returns acf array length N (float), with NaN where denom==0.
    This routine uses int64 accumulation to avoid overflow and then converts to float.
    """
    if arr.size == 0:
        return np.zeros(0, dtype=float)

    n_pairs, N = arr.shape
    acf = np.full(N, np.nan, dtype=float)

    # convert to int64 for safe sums
    arr64 = arr.astype(np.int64)

    # Precompute prefix sums per pair for fast denom? Not necessary here.
    # We'll compute numerator and denom for each tau.
    for tau in range(N):
        end = N - tau
        if end <= 0:
            acf[tau] = np.nan
            continue
        # slice arrays: left = arr64[:, :end], right = arr64[:, tau:tau+end]
        left = arr64[:, :end]
        right = arr64[:, tau:tau+end]
        # numerator: sum over p of dot(left[p], right[p])
        # vectorized: elementwise multiplication then sum
        prod = left * right  # shape (n_pairs, end)
        # sum using int64
        numer = int(prod.sum())   # safe int
        # denom: sum over p of left[p].sum()
        denom = int(left.sum())
        if denom == 0:
            acf[tau] = np.nan
        else:
            val = numer / denom
            # numerical safety: theoretically 0 <= val <= 1; enforce small tolerance
            if val < 0 and val > -1e-12:   # tiny negative due to fp
                val = 0.0
            # guard against small floating rounding slightly >1
            if val > 1 and val < 1 + 1e-12:
                val = 1.0
            # final clamp to [0,1]
            acf[tau] = min(max(val, 0.0), 1.0)
        if verbose and (tau % max(1, N//10) == 0):
            print(f"    intermittent ACF: tau {tau}/{N}")
    return acf


def analyze_hbonds_intermittent_blockavg(pkl_file, out_dir, skip_factor=1, nsteps_local=None, block_size=1000, verbose_local=True):
    with open(pkl_file, "rb") as f:
        hbonds_per_timestep_all = pickle.load(f)

    timesteps_all = sorted(hbonds_per_timestep_all.keys())
    timesteps = timesteps_all[::skip_factor]
    if nsteps_local is not None:
        timesteps = timesteps[:nsteps_local]
    effective_dt = dt * skip_factor
    N = len(timesteps)
    if N == 0:
        raise RuntimeError("No timesteps found!")

    nblocks = N // block_size
    if nblocks == 0:
        raise RuntimeError("Block size larger than trajectory length!")

    time = np.arange(block_size) * effective_dt

    results = {}
    all_cols = [time]
    header = ["time(ps)"]

    for bt in bond_types:
        if verbose_local:
            print(f"[INFO] Intermittent block-avg: {bt}, N={N}, blocks={nblocks}")

        # collect ACFs & taus for each block
        acf_blocks = []
        tau_blocks = []

        for b in range(nblocks):
            t_block = timesteps[b*block_size:(b+1)*block_size]
            sampled = {ts: hbonds_per_timestep_all[ts] for ts in t_block}
            arr, pairs_list = build_presence_arrays(sampled, t_block, bt)

            if arr.shape[0] == 0:
                continue

            acf = compute_intermittent_acf_from_arr(arr, verbose=False)
            if np.isfinite(acf[0]) and acf[0] != 0:
                acf = acf / acf[0]
            else:
                acf = np.nan_to_num(acf, nan=0.0)

            tau = float(np.trapz(np.nan_to_num(acf), time))
            acf_blocks.append(acf)
            tau_blocks.append(tau)

        if len(acf_blocks) == 0:
            acf_mean = np.full(block_size, np.nan)
            tau_mean, tau_std = np.nan, np.nan
        else:
            acf_mean = np.mean(acf_blocks, axis=0)
            tau_mean = float(np.mean(tau_blocks))
            tau_std = float(np.std(tau_blocks, ddof=1)) if len(tau_blocks) > 1 else 0.0

        results[bt] = {
            "acf_mean": acf_mean,
            "tau_mean": tau_mean,
            "tau_std": tau_std,
            "nblocks": nblocks
        }

        all_cols.append(acf_mean)
        header.append(f"{bt}_intermittent")

        # plot mean ACF
        plt.figure(figsize=(6,4))
        plt.plot(time, acf_mean, label=f"{bt} mean ACF", lw=1.5)
        plt.xlabel("Time (ps)"); plt.ylabel("C(t)")
        plt.title(f"{bt} intermittent ACF (block-avg)")
        plt.ylim(-0.05,1.05); plt.grid(alpha=0.2); plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"intermittent_blockavg_acf_{bt}.png"), dpi=300)
        plt.close()

    # save combined dat
    stacked = np.column_stack(all_cols)
    dat_file = os.path.join(out_dir, "hbonds_intermittent_blockavg_all.dat")
    np.savetxt(dat_file, stacked, header=" ".join(header), fmt="%.6g")

    # save summary
    summary_file = os.path.join(out_dir, "hbonds_intermittent_blockavg_summary.txt")
    with open(summary_file, "w") as fh:
        fh.write("# Intermittent ACF (block average) summary\n")
        fh.write(f"# dt={dt} ps, skip_factor={skip_factor}, effective_dt={effective_dt} ps\n")
        fh.write(f"# total timesteps used = {N}, block_size = {block_size}, nblocks = {nblocks}\n\n")
        for bt in bond_types:
            ent = results.get(bt, {})
            fh.write(f"{bt}:\n")
            fh.write(f"  tau_mean_ps: {ent.get('tau_mean')}\n")
            fh.write(f"  tau_std_ps: {ent.get('tau_std')}\n")
            fh.write(f"  nblocks: {ent.get('nblocks')}\n\n")

    if verbose_local:
        print(f"[INFO] Block-avg results saved to {summary_file}")
    return results


if __name__ == "__main__":
    res_dir = r"C:/LEMTA/Data/SiO2/ClayFF/Hydrogen_Network/bond_lifetime/100_perc/Results/"
    pkl_file = os.path.join(res_dir, "hbonds_per_timestep.pkl")

    results = analyze_hbonds_intermittent_blockavg(
        pkl_file, res_dir,
        skip_factor=skip_factor,
        nsteps_local=nsteps,
        block_size=block_size,
        verbose_local=True
    )