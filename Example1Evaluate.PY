"""
validate_CL_metrics_1.py
========================
Python translation of validate_evaluate_CL_metrics.m (Validation Script 1).

Generates a piecewise synthetic MAP signal with analytically known
metric values and verifies that evaluate_CL_metrics produces matching
results.

Signal design (target = 70 mmHg, 100 Hz, 60 min):
  Phase 1  t = 0 – 5 min   Linear rise  42 → 84 mmHg
  Phase 2  t = 5 – 10 min  Linear fall  84 → 70 mmHg
  Phase 3  t = 10 – 60 min Constant     70 mmHg

USAGE
-----
    python validate_CL_metrics_1.py

Place evaluate_CL_metrics.py in the same directory.

REQUIREMENTS
------------
    pip install numpy pandas scipy matplotlib
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from evaluate_CL_metrics import evaluate_CL_metrics

# -----------------------------------------------------------------------
# 1. INPUTS
# -----------------------------------------------------------------------
target       = 70.0
Dt           = 0.01                          # min (100 Hz)
Times        = np.arange(0, 60 + Dt/2, Dt)  # 6001 samples
N            = len(Times)                    # 6001
ind_CL_start = 1                             # 1-based (t = 0 min)
ind_CL_end   = N                             # 1-based (t = 60 min)
Hemorrhage   = np.zeros(N)                   # no disturbance → SSE Case 1

# -----------------------------------------------------------------------
# 2. SIGNAL CONSTRUCTION
# -----------------------------------------------------------------------
# MATLAB:  MAP(1:501)    = linspace(42, 84, 501)   → Python: MAP[0:501]
# MATLAB:  MAP(501:1001) = linspace(84, 70, 501)   → Python: MAP[500:1001]
# MATLAB:  MAP(1001:N)   = 70                      → Python: MAP[1000:]
MAP          = np.zeros(N)
MAP[0:501]   = np.linspace(42, 84, 501)   # Phase 1: rise
MAP[500:1001] = np.linspace(84, 70, 501)  # Phase 2: fall  (MAP[500]=84, same)
MAP[1000:]   = 70.0                        # Phase 3: steady state

# -----------------------------------------------------------------------
# 3. INDEPENDENT EXPECTED VALUE DERIVATIONS
# -----------------------------------------------------------------------
s   = ind_CL_start - 1    # 0-based start
e   = ind_CL_end          # Python exclusive end
PE_ind  = (MAP[s:e] - target) / target * 100.0
n_eng   = e - s
MAP_eng = MAP[s:e]

# (i) MDPE  — Phase 3 has 5001/6001 samples at PE=0 → median = 0%
MDPE_exp  = float(np.median(PE_ind))          # 0.00 %

# (ii) MDAPE — same reasoning → 0%
MDAPE_exp = float(np.median(np.abs(PE_ind)))  # 0.00 %

# (iii) Wobble = median(|PE - MDPE|) = median(|PE|) = 0%
Wobble_exp = float(np.median(np.abs(PE_ind - MDPE_exp)))

# (iv) POT
per10_exp = float(np.sum(np.abs(MAP_eng - target) < 0.1*target) * 100 / n_eng)
per20_exp = float(np.sum(np.abs(MAP_eng - target) < 0.2*target) * 100 / n_eng)
per30_exp = float(np.sum(np.abs(MAP_eng - target) < 0.3*target) * 100 / n_eng)

# (v) Divergence  (OLS slope)
n_pe   = len(PE_ind)
t_temp = np.arange(n_pe) * Dt
mean_t = float(np.mean(t_temp))
divergence_exp = float(
    (np.sum(np.abs(PE_ind) * t_temp) - np.sum(np.abs(PE_ind)) * mean_t)
    / (np.sum(t_temp**2) - np.sum(t_temp)**2 / n_pe)
)

# (vi) Rise time
# MAP_focus[0]=42, T10=44.8, T90=67.2, Phase 1 step=0.084
# idx_T10=34  (MAP[33]=44.772 → round=44.8 ≥ 44.8 ✓)
# idx_T90=301 (MAP[300]=67.2  → round=67.2 ≥ 67.2 ✓)
# risetime = (301-34)*0.01 = 2.67 min
MAP_focus = MAP[s:e]
T10_lvl   = round(MAP_focus[0] + 0.1*(target - MAP_focus[0]), 1)   # 44.8
T90_lvl   = round(MAP_focus[0] + 0.9*(target - MAP_focus[0]), 1)   # 67.2
r_map     = np.round(MAP_focus, 1)
idx_T10_exp = int(np.where(r_map >= T10_lvl)[0][0])   # 0-based → 33
idx_T90_exp = int(np.where(r_map >= T90_lvl)[0][0])   # 0-based → 300
risetime_exp = float((idx_T90_exp - idx_T10_exp) * Dt)  # 2.67 min

# (vii) Overshoot  max=84 → (84-70)/70*100 = 20%
overshoot_exp = float((np.max(MAP_focus) - target) / target * 100)

# (viii) Settling time
# Settling uses ≤ 0.1 (inclusive), per10 uses strict < 0.1.
# MAP(750) = 84-250*0.028=77.0  PF=7/70=0.1 ≤ 0.1 ✓ → settling starts here
# MATLAB Times(751) = 750*0.01 = 7.50 min  →  Python Times[750] = 7.50 min
settle_exp = float(Times[750])   # 7.50 min

# (ix) SSE  — no hemorrhage → Case 1 → last 2 min at MAP=70 → SSE=0
w         = round(2.0 / Dt)      # 200 samples
SSE_exp   = float(np.mean(MAP[ind_CL_end - w : ind_CL_end]) - target)  # 0.00

# -----------------------------------------------------------------------
# 4. RUN evaluate_CL_metrics
# -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end,
                              MAP, Times, target, Hemorrhage)

# -----------------------------------------------------------------------
# 5. COMPARISON TABLE
# -----------------------------------------------------------------------
tol = 1e-6

print()
print('=' * 60)
print('  Validation Script 1  —  Results')
print('=' * 60)
print(f"  {'Metric':<22}  {'Expected':>10}  {'Computed':>10}  {'Status'}")
print('-' * 60)

labels = ['MDPE (%)', 'MDAPE (%)', 'Wobble (%)',
          'per10 (%)', 'per20 (%)', 'per30 (%)',
          'risetime (min)', 'overshoot (%)', 'settlingtime (min)',
          'divergence (%/min)', 'SSE (mmHg)']

exp_vals  = [MDPE_exp, MDAPE_exp, Wobble_exp,
             per10_exp, per20_exp, per30_exp,
             risetime_exp, overshoot_exp, settle_exp,
             divergence_exp, SSE_exp]

keys = ['MDPE', 'MDAPE', 'Wobble',
        'per10', 'per20', 'per30',
        'risetime', 'overshoot', 'settlingtime',
        'divergence', 'SSE']

all_pass = True
for label, exp, key in zip(labels, exp_vals, keys):
    comp = result[key]
    ok   = (not math.isnan(exp) and not math.isnan(comp)
            and abs(exp - comp) < tol)
    status = 'PASS' if ok else f'FAIL (diff={exp-comp:.2e})'
    if not ok:
        all_pass = False
    print(f"  {label:<22}  {exp:>10.4f}  {comp:>10.4f}  {status}")

print('-' * 60)
if all_pass:
    print(f'  All metrics matched within tolerance ({tol:.0e}).')
else:
    print('  One or more metrics did NOT match.')
print('=' * 60)

# -----------------------------------------------------------------------
# 6. PLOT
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(Times, MAP, 'b-', linewidth=1.5, label='MAP signal')
ax.axhline(y=target,       color='k', linestyle='--', linewidth=1.2,
           label='Target (70 mmHg)')
ax.axhline(y=target*1.1,   color='r', linestyle=':',  linewidth=1.0)
ax.axhline(y=target*0.9,   color='r', linestyle=':',  linewidth=1.0,
           label='±10% band')
ax.axvline(x=Times[idx_T10_exp], color='m', linestyle='--', linewidth=1.0)
ax.axvline(x=Times[idx_T90_exp], color='m', linestyle='--', linewidth=1.0)
ax.axvline(x=result['settlingtime'], color='g', linestyle='--',
           linewidth=1.2,
           label=f"Settle = {result['settlingtime']:.2f} min")
ax.text(Times[idx_T10_exp]+0.2, T10_lvl+2, 'T₁₀', color='m', fontsize=9)
ax.text(Times[idx_T90_exp]+0.2, T90_lvl+2, 'T₉₀', color='m', fontsize=9)
ax.text(30, 82, f"Overshoot = {result['overshoot']:.1f}%",
        color='r', fontsize=9)
ax.set_xlabel('Time (min)')
ax.set_ylabel('MAP (mmHg)')
ax.set_title('Validation Script 1 — Synthetic Signal')
ax.legend(loc='lower right')
ax.set_ylim([35, 92])
ax.grid(True)
plt.tight_layout()
plt.show()
