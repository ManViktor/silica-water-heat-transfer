# -*- coding: utf-8 -*-
"""
analyze_hbonds_continuous_block.py
Continuous ACF for hydrogen bonds, block-averaging version.
- Splits trajectory into blocks (block_size steps)
- Computes ACF for each block in parallel
- Normalizes, integrates/fits, averages over blocks
- Saves .dat (average ACFs), summary (mean+std tau), and per-type plots

Author: adapted for Viktor
"""

import os
import pickle
import numpy as np
from scipy.optimize import curve_fit
import multiprocessing as mp
import gc  # for garbage collection

# Set matplotlib backend to Agg for headless environments (no DISPLAY)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------- USER SETTINGS ----------------
dt = 0.01           # ps between saved frames
skip_factor = 1
nsteps = 10000       # total steps to use (reduced for memory test; increase if OK)
block_size = 200    # steps per block (reduced for memory)
bond_types = ["OH_donor", "OH_acceptor", "CH3_donor"]
max_fit_time = None
verbose = True
num_processes = 1   # number of CPU cores to use
# ------------------------------------------------

def exp_decay(t, tau):
    return np.exp(-t / tau)

def build_presence_arrays(hbonds_per_timestep, timesteps, bt):
    pairs_set = set()
    for ts in timesteps:
        entry = hbonds_per_timestep[ts]
        current = entry.get(bt, set()) if isinstance(entry, dict) else (entry if bt == "all" else set())
        pairs_set.update(current)
    pairs_list = sorted(pairs_set)
    arr = np.zeros((len(pairs_list), len(timesteps)), dtype=np.int8)
    for j, ts in enumerate(timesteps):
        entry = hbonds_per_timestep[ts]
        current = entry.get(bt, set()) if isinstance(entry, dict) else (entry if bt == "all" else set())
        for i, hb in enumerate(pairs_list):
            if hb in current:
                arr[i, j] = 1
    return {pair: arr[i] for i, pair in enumerate(pairs_list)}

def compute_continuous_acf_from_presence(presence_dict, verbose=False):
    pairs = list(presence_dict.keys())
    if len(pairs) == 0: 
        return np.zeros(0)
    n = len(next(iter(presence_dict.values())))
    acf = np.full(n, np.nan, dtype=float)
    for tau in range(n):
        numer, denom = 0, 0
        end = n - tau
        if end <= 0: 
            continue
        window_len = tau + 1
        for p in pairs:
            arr = presence_dict[p]
            denom += arr[:end].sum()
            if window_len == 1:
                numer += arr[:end].sum()
            else:
                win = np.convolve(arr, np.ones(window_len, dtype=int), mode='full')[:n]
                numer += np.sum(win[:end] == window_len)
        acf[tau] = (numer / denom) if denom > 0 else np.nan
        if verbose and tau % max(1, n//10) == 0:
            print(f"    continuous ACF: tau={tau}/{n}")
    return acf

def safe_exponential_fit(time, acf, max_time=None):
    mask = np.isfinite(acf) & (acf > 0)
    if not np.any(mask):
        return np.nan, {"method":"none"}
    max_idx = np.where(acf >= 0.01)[0][-1] if np.any(acf >= 0.01) else 3
    if max_time is not None and len(time) > 1:
        max_idx_time = int(max_time / (time[1]-time[0]))
        max_idx = min(max_idx, max_idx_time)
    x, y = time[:max_idx+1], acf[:max_idx+1]
    tau0 = max(1e-3, np.trapz(y, x) / (y[0] if y[0]>0 else 1.0))
    try:
        popt,_ = curve_fit(exp_decay, x, y, p0=[tau0], maxfev=10000)
        tau_fit = float(popt[0])
        return tau_fit, {"method":"exp_fit"}
    except Exception as e:
        return np.trapz(y, x), {"method":f"integral_fallback ({str(e)})"}

def compute_block_acf(args):
    """Compute ACF for one block. Args: (hbonds_all, tblock, bt, block_size, effective_dt, max_fit_time)"""
    hbonds_all, tblock, bt, block_size, effective_dt, max_fit_time = args
    sampled = {ts: hbonds_all[ts] for ts in tblock}
    presence = build_presence_arrays(sampled, tblock, bt)
    acf = compute_continuous_acf_from_presence(presence, verbose=False)
    if acf.size == 0: 
        return None, np.nan
    time = np.arange(block_size) * effective_dt
    acf = acf / acf[0] if acf[0] != 0 else np.nan_to_num(acf)
    tau, _ = safe_exponential_fit(time, acf, max_time=max_fit_time)
    gc.collect()  # Clean up memory after block
    return acf, tau

def analyze_continuous_block(pkl_file, out_dir):
    with open(pkl_file,"rb") as f:
        hbonds_per_timestep_all = pickle.load(f)
    timesteps_all = sorted(hbonds_per_timestep_all.keys())
    timesteps = timesteps_all[::skip_factor][:nsteps]
    effective_dt = dt * skip_factor
    N = len(timesteps)
    if N == 0:
        raise RuntimeError("No timesteps found!")
    time = np.arange(block_size) * effective_dt

    results = {}
    all_cols = [time]
    header = ["time(ps)"]

    nblocks = N // block_size
    print(f"[INFO] Splitting into {nblocks} blocks of {block_size} steps")
    print(f"[INFO] Using {num_processes} processes for parallel computation")

    # Prepare pool
    with mp.Pool(processes=num_processes) as pool:
        for bt in bond_types:
            print(f"[INFO] Processing {bt}...")
            acfs_blocks = []
            taus = []
            
            # Prepare args for each block
            block_args = [
                (hbonds_per_timestep_all, timesteps[b*block_size:(b+1)*block_size], bt, block_size, effective_dt, max_fit_time)
                for b in range(nblocks)
            ]
            
            # Parallel compute
            results_blocks = pool.map(compute_block_acf, block_args)
            
            for acf, tau in results_blocks:
                if acf is not None:
                    acfs_blocks.append(acf)
                    taus.append(tau)
            
            if len(acfs_blocks) == 0:
                mean_acf = np.full(block_size, np.nan)
                tau_mean = np.nan
                tau_std = np.nan
            else:
                mean_acf = np.nanmean(acfs_blocks, axis=0)
                tau_mean = float(np.nanmean(taus))
                tau_std = float(np.nanstd(taus))
            results[bt] = {"acf": mean_acf, "tau_mean": tau_mean, "tau_std": tau_std}

            all_cols.append(mean_acf)
            header.append(f"{bt}_continuous")

            # plot per type
            plt.figure(figsize=(6,4))
            plt.plot(time, mean_acf, label=f"{bt} mean ACF")
            if np.isfinite(tau_mean):
                plt.plot(time, exp_decay(time,tau_mean),'--',label=f"fit t={tau_mean:.3f} ps")
            plt.xlabel("Time (ps)"); plt.ylabel("C(t)"); plt.title(f"{bt} continuous ACF (block avg)")
            plt.legend(); plt.tight_layout()
            plt.savefig(os.path.join(out_dir,f"continuous_acf_{bt}.png"),dpi=300); plt.close()

            print(f"  {bt}: tau_mean={tau_mean:.4g} ps, std={tau_std:.4g} ps")

    # save combined .dat
    stacked = np.column_stack(all_cols)
    dat_file = os.path.join(out_dir,"hbonds_acf_continuous_all.dat")
    np.savetxt(dat_file, stacked, header=" ".join(header), fmt="%.6g")
    print(f"[INFO] All continuous ACFs saved to {dat_file}")

    # save summary
    summary_file = os.path.join(out_dir,"hbonds_acf_continuous_summary.txt")
    with open(summary_file,"w") as fh:
        fh.write("# Continuous ACF summary (block averaged)\n")
        fh.write(f"# dt={dt} ps, skip_factor={skip_factor}, effective_dt={effective_dt} ps\n")
        fh.write(f"# total_timesteps_used={N}, block_size={block_size}, nblocks={nblocks}\n")
        fh.write(f"# num_processes={num_processes}\n\n")
        for bt in bond_types:
            ent = results[bt]
            fh.write(f"{bt}:\n")
            fh.write(f"  tau_mean_ps: {ent['tau_mean']}\n")
            fh.write(f"  tau_std_ps: {ent['tau_std']}\n\n")

    # plot combined
    plt.figure(figsize=(8,5))
    for bt in bond_types:
        acf = results[bt]["acf"]
        if acf is not None and acf.size>0:
            plt.plot(time, acf, label=f"{bt}")
    plt.xlabel("Time (ps)"); plt.ylabel("C(t)"); plt.title("Continuous ACF (all bond types, block avg)")
    plt.ylim(-0.05,1.05); plt.legend(); plt.tight_layout()
    png_all = os.path.join(out_dir,"continuous_acf_all.png")
    plt.savefig(png_all,dpi=300); plt.close()
    print(f"[INFO] Combined plot saved to {png_all}")

    gc.collect()  # Final cleanup
    return results

if __name__=="__main__":
    # Set multiprocessing start method for compatibility (esp. on clusters)
    mp.set_start_method('spawn', force=True)
    pkl_file = "hbonds_per_timestep.pkl"
    out_dir = "."
    analyze_continuous_block(pkl_file, out_dir)