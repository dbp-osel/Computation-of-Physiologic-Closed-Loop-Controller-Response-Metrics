"""
validate_CL_metrics_2.py
========================
Python translation of validate_evaluate_CL_metrics_2.m (Validation Script 2).

Challenge signal exercising three code paths not covered by Scripts 1–2:
  (A) ind_CL_start ≠ 1  — pre-engagement period at MAP=45 mmHg
  (B) SSE Case 2        — hemorrhage within 10 min of CL end
  (C) Three settling groups — two fail, third qualifies

Signal design (target = 70 mmHg, 100 Hz, 70 min):
  Pre-CL   t =  0–10 min   Constant       45 mmHg  (not engaged)
  Phase 1  t = 10–14 min   Linear rise    45 → 81 mmHg
  Phase 2  t = 14–18 min   Linear fall    81 → 61 mmHg
  Phase 3  t = 18–22 min   Linear rise    61 → 73 mmHg
  Phase 4  t = 22–55 min   Constant       73 mmHg
  Phase 5  t = 55–63 min   Hemorrhage     73 → 55 mmHg
  Phase 6  t = 63–70 min   Recovery       55 → 70 mmHg

USAGE
-----
    python validate_CL_metrics_2.py

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
Dt           = 0.01
Times        = np.arange(0, 70 + Dt/2, Dt)   # 7001 samples
N            = len(Times)                     # 7001
ind_CL_start = 1001                           # 1-based (t = 10 min)
ind_CL_end   = N                              # 1-based (t = 70 min)

# -----------------------------------------------------------------------
# 2. SIGNAL CONSTRUCTION
# -----------------------------------------------------------------------
# Index mapping:  MATLAB MAP(k) [1-based]  →  Python MAP[k-1] [0-based]
# MATLAB MAP(1:1001)    = 45               → MAP[0:1001]
# MATLAB MAP(1001:1401) = linspace(...)    → MAP[1000:1401]
# etc.
MAP = np.zeros(N)
MAP[0:1001]    = 45.0                          # Pre-CL (t=0-10 min)
MAP[1000:1401] = np.linspace(45, 81, 401)     # Phase 1: rise  step=0.09
MAP[1400:1801] = np.linspace(81, 61, 401)     # Phase 2: fall  step=0.05
MAP[1800:2201] = np.linspace(61, 73, 401)     # Phase 3: recovery
MAP[2200:5501] = 73.0                          # Phase 4: steady state
MAP[5500:6301] = np.linspace(73, 55, 801)     # Phase 5: hemorrhage drop
MAP[6300:N]    = np.linspace(55, 70, N-6300)  # Phase 6: recovery

Hemorrhage = np.zeros(N)
Hemorrhage[5500:6301] = -0.1    # hemorrhage active t=55-63 min

# -----------------------------------------------------------------------
# 3. INDEPENDENT EXPECTED VALUE DERIVATIONS
# -----------------------------------------------------------------------
s       = ind_CL_start - 1    # 0-based: 1000
e       = ind_CL_end          # Python exclusive end: 7001
PE_ind  = (MAP[s:e] - target) / target * 100.0
n_eng   = e - s               # 6001
MAP_eng = MAP[s:e]

# (i) MDPE — Phase 4 has 3301/6001 = 55% at PE=4.286% → median = 4.286%
MDPE_exp  = float(np.median(PE_ind))

# (ii) MDAPE
MDAPE_exp = float(np.median(np.abs(PE_ind)))

# (iii) Wobble — Phase 4 constant 73 gives |PE-MDPE|=0 for 55% → Wobble=0%
Wobble_exp = float(np.median(np.abs(PE_ind - MDPE_exp)))

# (iv) POT
per10_exp = float(np.sum(np.abs(MAP_eng-target) < 0.1*target)*100/n_eng)
per20_exp = float(np.sum(np.abs(MAP_eng-target) < 0.2*target)*100/n_eng)
per30_exp = float(np.sum(np.abs(MAP_eng-target) < 0.3*target)*100/n_eng)

# (v) Divergence
n_pe   = len(PE_ind)
t_temp = np.arange(n_pe) * Dt
mean_t = float(np.mean(t_temp))
divergence_exp = float(
    (np.sum(np.abs(PE_ind)*t_temp) - np.sum(np.abs(PE_ind))*mean_t)
    / (np.sum(t_temp**2) - np.sum(t_temp)**2/n_pe)
)

# (vi) Rise time
# MAP_focus[0] = MAP[1000] = 45 mmHg  (ind_CL_start-1 = 1000)
# T10=47.5, T90=67.5, Phase 1 step=0.09
# MATLAB idx_T10=29 (0-based: 28)  MAP[1028]=45+28*0.09=47.52 → round=47.5 ✓
# MATLAB idx_T90=251(0-based: 250) MAP[1250]=45+250*0.09=67.5 → round=67.5 ✓
# risetime = (250-28)*0.01 = 2.22 min
MAP_focus = MAP[s:e]    # local engagement array (0-based locally)
T10_lvl   = round(MAP_focus[0] + 0.1*(target - MAP_focus[0]), 1)   # 47.5
T90_lvl   = round(MAP_focus[0] + 0.9*(target - MAP_focus[0]), 1)   # 67.5
r_map     = np.round(MAP_focus, 1)
idx_T10_exp = int(np.where(r_map >= T10_lvl)[0][0])   # local 0-based = 28
idx_T90_exp = int(np.where(r_map >= T90_lvl)[0][0])   # local 0-based = 250
risetime_exp = float((idx_T90_exp - idx_T10_exp) * Dt)  # 2.22 min

# (vii) Overshoot  max=81 → (81-70)/70*100 = 15.714%
overshoot_exp = float((np.max(MAP_focus) - target) / target * 100)

# (viii) Settling time — three groups
#
# Group 1 (Phase 1 upswing, local idx 210–336):
#   MAP[1210]=45+210*0.09=63.9... let me derive:
#   MAP enters band: 45+(k)*0.09>=63 → k>=200 → local idx 200 (0-based)
#   MAP[200]=45+200*0.09=63.0  PF=7/70=0.1 ≤ 0.1 ✓
#   MAP exits band: (k)*0.09>32 → k>355.6 → local idx 356
#   MAP[356]=45+356*0.09=77.04  PF>0.1 ✗  last in-band: local 355
#   Duration=(355-200)*0.01=1.55 min < 10 → FAIL
#
# Group 2 (Phase 2 downswing, local 480–760):
#   Phase 2: MAP=81-(k-400)*0.05, enters ≤77 at k=400+80=480
#   MAP[480]=81-80*0.05=77.0  PF=0.1 ✓
#   Exits below 63: (k-400)*0.05>18 → k>760  MAP[760]=63.0 ✓
#   MAP[761]=62.95  PF>0.1 ✗  last in-band: local 760
#   Duration=(760-480)*0.01=2.80 min < 10 → FAIL
#
# Group 3 (Phase 3 recovery + Phase 4 + into Phase 5, qualifies):
#   Phase 3 enters band at local idx 867:
#   61+(k-800)*0.03>=63 → k>=800+66.7 → k=867
#   MAP[867]=61+67*0.03=63.01  PF≤0.1 ✓  MAP[866]=62.98  PF>0.1 ✗
#   Phase 5 drops below 63 at:
#   73-(k-4500)*0.0225<63 → k>4500+444.4 → k=4945
#   MAP[4944]=73-444*0.0225=63.01  PF≤0.1 ✓  MAP[4945]=62.99  PF>0.1 ✗
#   Duration=(4944-867)*0.01=40.77 min ≥ 10 → QUALIFIES
#   settle = t_local[867] = Times[s+867] = Times[1000+867] = Times[1867]
#   Times[1867] = 1867*0.01 = 18.67 min
settle_exp = float(Times[1867])    # 18.67 min

# (ix) SSE — Case 2
# Hem ends at global idx ~6310 (after smoothing).
# (7001-6310)*0.01 ≈ 6.9 min < 10 → Case 2.
# While loop converges before hemorrhage onset (~global 5500).
# 2-min window: global 5301-5500 → Phase 4 (MAP=73) → SSE=3.00
SSE_exp = 3.00

# -----------------------------------------------------------------------
# 4. RUN evaluate_CL_metrics
# -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end,
                              MAP, Times, target, Hemorrhage)

# -----------------------------------------------------------------------
# 5. COMPARISON TABLE
# -----------------------------------------------------------------------
tol = 1e-4    # slightly relaxed for SSE Case 2 smoothing variation

print()
print('=' * 65)
print('  Validation Script 2  —  Challenge Signal')
print(f'  ind_CL_start={ind_CL_start} (t=10 min) | SSE Case 2 | 3 groups')
print('=' * 65)
print(f"  {'Metric':<22}  {'Expected':>10}  {'Computed':>10}  {'Status'}")
print('-' * 65)

labels = ['MDPE (%)', 'MDAPE (%)', 'Wobble (%)',
          'per10 (%)', 'per20 (%)', 'per30 (%)',
          'risetime (min)', 'overshoot (%)', 'settlingtime (min)',
          'divergence (%/min)', 'SSE (mmHg)']

exp_vals = [MDPE_exp, MDAPE_exp, Wobble_exp,
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

print('-' * 65)
if all_pass:
    print(f'  All metrics matched within tolerance ({tol:.0e}).')
else:
    print('  One or more metrics did NOT match.')
print('=' * 65)

# --- Challenge A: show what wrong PE indexing would produce ---
print()
print('--- Challenge A: PE indexing ---')
PE_buggy   = (MAP[0:n_eng] - target) / target * 100.0
MDPE_buggy = float(np.median(PE_buggy))
print(f'  Correct MDPE (MAP[1000:7001]): {result["MDPE"]:7.3f} %')
print(f'  Buggy   MDPE (MAP[0:6001]):    {MDPE_buggy:7.3f} %  ← pre-CL dominates')

# --- Challenge B: SSE Case 2 ---
print()
print('--- Challenge B: SSE Case 2 ---')
print(f'  Hemorrhage ends at t≈63 min; CL ends at t=70 min')
print(f'  Post-hemorrhage window ≈ 7 min < 10 min → Case 2 triggered')
print(f'  SSE = {result["SSE"]:.2f} mmHg  (window lands in Phase 4, MAP=73)')

# --- Challenge C: settling groups ---
print()
print('--- Challenge C: settling groups ---')
print('  Group 1 (Phase 1):  ~1.55 min  < 10  → FAIL')
print('  Group 2 (Phase 2):  ~2.80 min  < 10  → FAIL')
print(f'  Group 3 (Phase 3+): ~40.77 min ≥ 10  → PASS → settle={result["settlingtime"]:.2f} min')

# -----------------------------------------------------------------------
# 6. PLOT
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5))

# Shade regions
ax.axvspan(0, 10, alpha=0.08, color='gray',  label='Pre-CL (excluded)')
ax.axvspan(10, 70, alpha=0.04, color='blue')
ax.axvspan(55, 63, alpha=0.15, color='red',  label='Hemorrhage (t=55–63)')

ax.plot(Times, MAP, 'b-', linewidth=1.5, label='MAP signal')
ax.axhline(y=target,     color='k', linestyle='--', linewidth=1.2,
           label='Target (70 mmHg)')
ax.axhline(y=target*1.1, color='r', linestyle=':',  linewidth=1.0)
ax.axhline(y=target*0.9, color='r', linestyle=':',  linewidth=1.0,
           label='±10% band')
ax.axvline(x=10,  color='k', linestyle='-.', linewidth=1.0,
           label='CL start (t=10)')
ax.axvline(x=Times[s + idx_T10_exp], color='m', linestyle='--', linewidth=1.0)
ax.axvline(x=Times[s + idx_T90_exp], color='m', linestyle='--', linewidth=1.0)
ax.axvline(x=result['settlingtime'], color='g', linestyle='--',
           linewidth=1.2,
           label=f"Settle = {result['settlingtime']:.2f} min (Group 3)")

ax.text(Times[s+idx_T10_exp]+0.2, T10_lvl+4, 'T₁₀', color='m', fontsize=8)
ax.text(Times[s+idx_T90_exp]+0.2, T90_lvl+4, 'T₉₀', color='m', fontsize=8)
ax.text(12, 88, f"Overshoot = {result['overshoot']:.2f}%",
        color='r', fontsize=8)
ax.text(12, 84, f"MDPE = {result['MDPE']:.3f}%", color='b', fontsize=8)
ax.text(12, 80, f"SSE  = {result['SSE']:.2f} mmHg (Case 2)",
        color='b', fontsize=8)
ax.text(3.5, 52, 'Pre-CL', color='gray', fontsize=8, fontstyle='italic')

ax.set_xlabel('Time (min)')
ax.set_ylabel('MAP (mmHg)')
ax.set_title('Validation Script 2 — Challenge Signal')
ax.legend(loc='lower right', fontsize=8)
ax.set_ylim([38, 95])
ax.grid(True)
plt.tight_layout()
plt.show()
