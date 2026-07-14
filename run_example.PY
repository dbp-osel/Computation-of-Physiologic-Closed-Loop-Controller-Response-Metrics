"""
run_example.py
==============
Python equivalent of run_example.m.

Loads test_data.mat, calls evaluate_CL_metrics, prints results,
and reproduces the two-panel figure (MAP vs Time / Hemorrhage vs Time).

USAGE
-----
    python run_example.py

Place test_data.mat, evaluate_CL_metrics.py, and this file in the
same directory before running.

REQUIREMENTS
------------
    pip install numpy pandas scipy matplotlib
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# ---- Load the .mat file -----------------------------------------------
# scipy.io.loadmat handles MATLAB files up to v7.2.
# If the file is MATLAB 7.3 (HDF5 format), h5py is used as a fallback.
try:
    import scipy.io
    data = scipy.io.loadmat('test_data.mat')

except Exception:
    # MATLAB 7.3 / HDF5 fallback
    import h5py
    data = {}
    with h5py.File('test_data.mat', 'r') as f:
        for k in f.keys():
            data[k] = np.array(f[k])

# ---- Extract variables ------------------------------------------------
# MATLAB saves scalars as 2-D arrays, so .flatten()[0] is used everywhere.
ind_CL_start      = int(   np.asarray(data['ind_CL_start']).flatten()[0])
ind_CL_end        = int(   np.asarray(data['ind_CL_end']).flatten()[0])
MAP_Values        = np.asarray(data['MAP_Values']).flatten().astype(float)
MAP_Times         = np.asarray(data['MAP_Times']).flatten().astype(float)
target            = float( np.asarray(data['target']).flatten()[0])
Hemorrhage_Values = np.asarray(data['Hemorrhage_Values']).flatten().astype(float)

# ---- Run the function -------------------------------------------------
from evaluate_CL_metrics import evaluate_CL_metrics

PCLC_char = evaluate_CL_metrics(
    ind_CL_start,
    ind_CL_end,
    MAP_Values,
    MAP_Times,
    target,
    Hemorrhage_Values
)

# ---- Print results (matches MATLAB's unsilenced struct display) --------
print('\nPCLC_char =')
print(f"  {'Metric':<15}  {'Value'}")
print(f"  {'-'*15}  {'-'*12}")
for key, val in PCLC_char.items():
    if isinstance(val, float) and math.isnan(val):
        print(f"  {key:<15}  NaN")
    else:
        print(f"  {key:<15}  {val:.4f}")
print()

# ---- Figure -----------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
fig.subplots_adjust(hspace=0.08)

# --- Top panel: MAP vs Time ---
ax1 = axes[0]
ax1.plot(MAP_Times, MAP_Values, 'b', linewidth=1.5, label='MAP')
ax1.axhline(y=target, color='r', linestyle='--', linewidth=2, label='Target')
ax1.axvline(
    x=MAP_Times[ind_CL_start - 1],   # 1-based → 0-based index
    color='g', linestyle='--', linewidth=1.5, label='PCLC Start'
)
ax1.axvline(
    x=MAP_Times[ind_CL_end - 1],
    color='m', linestyle='--', linewidth=1.5, label='PCLC End'
)
ax1.set_ylabel('MAP (mmHg)')
ax1.legend(loc='best')
ax1.grid(True)

# --- Bottom panel: Hemorrhage vs Time ---
# MATLAB stairs → matplotlib step with where='post'
ax2 = axes[1]
ax2.step(MAP_Times, Hemorrhage_Values, 'k', linewidth=1.5,
         where='post')
ax2.axvline(
    x=MAP_Times[ind_CL_start - 1],
    color='g', linestyle='--', linewidth=1.5, label='PCLC Start'
)
ax2.axvline(
    x=MAP_Times[ind_CL_end - 1],
    color='m', linestyle='--', linewidth=1.5, label='PCLC End'
)
ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Hemorrhage (l/min)')
ax2.grid(True)

plt.tight_layout()
plt.show()
