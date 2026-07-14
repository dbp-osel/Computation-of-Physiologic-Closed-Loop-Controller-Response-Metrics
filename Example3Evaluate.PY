"""
validate_CL_metrics_3.py
========================
Python translation of validate_evaluate_CL_metrics_4.m (Validation Script 3).

Maximum challenge signal — 3 hemorrhages, 4 settling groups, Wobble≠0.

  (A) ind_CL_start = 1001  — pre-engagement at MAP=42 mmHg
  (B) THREE hemorrhage events:
        Hem1  t=24-26 min  (mid-engagement)
        Hem2  t=42-45 min  (mid-engagement)
        Hem3  t=88-92 min  (near end → SSE Case 2)
  (C) FOUR settling groups:
        Group 1   1.34 min  FAIL  (trivially short)
        Group 2   8.70 min  FAIL  (only 1.30 min short of threshold!)
        Group 3  15.62 min  PASS  → settle = 28.25 min
        Group 4  42.97 min  PASS  (subsequent)
  (D) Non-zero Wobble — sinusoidal Phase 7 (73±4 mmHg, period=8 min)

Signal design (target=70 mmHg, 100 Hz, 100 min):
  Pre-CL    t =  0-10 min   Constant   42 mmHg
  Phase 1   t = 10-14 min   Rise       42 → 84 mmHg
  Phase 2   t = 14-18 min   Fall       84 → 73 mmHg
  Phase 3   t = 18-24 min   Constant   73 mmHg
  Hem1      t = 24-26 min   Fall       73 → 57 mmHg
  Phase 4   t = 26-32 min   Rise       57 → 73 mmHg
  Phase 5   t = 32-42 min   Constant   73 mmHg
  Hem2      t = 42-45 min   Fall       73 → 57 mmHg
  Phase 6   t = 45-51 min   Rise       57 → 73 mmHg
  Phase 7   t = 51-67 min   Sinusoidal 73+4*sin(2π(t-51)/8)
  Phase 8   t = 67-85 min   Constant   73 mmHg
  Phase 8.5 t = 85-88 min   Constant   73 mmHg  (SSE window)
  Hem3      t = 88-92 min   Fall       73 → 55 mmHg
  Phase 9   t = 92-100 min  Rise       55 → 70 mmHg

Analytically known values:
  MDPE        =  4.286 %    (73-70)/70*100
  Wobble      >  0 %        sinusoidal phase, <42% samples at exactly 73
  Overshoot   = 20.000 %    (84-70)/70*100
  Rise time   =  2.13 min   (idx_T90-idx_T10)*0.01 = (241-28)*0.01
  Settle time = 28.25 min   Times[2825]  (Group 3, local idx 1826)
  SSE         =  3.000 mmHg Case 2, window before Hem1 at Phase 3 MAP=73

USAGE
-----
    python validate_CL_metrics_3.py

REQUIREMENTS
------------
    pip install numpy pandas scipy matplotlib
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from evaluate_CL_metrics import evaluate_CL_metrics

# -----------------------------------------------------------------------
# 1. INPUTS
# -----------------------------------------------------------------------
target       = 70.0
Dt           = 0.01
Times        = np.arange(0, 100 + Dt/2, Dt)  # 10001 samples
N            = len(Times)                     # 10001
ind_CL_start = 1001                           # 1-based (t = 10 min)
ind_CL_end   = N                              # 1-based (t = 100 min)

# -----------------------------------------------------------------------
# 2. SIGNAL CONSTRUCTION
# -----------------------------------------------------------------------
# Index rule:  MATLAB MAP(a:b) [1-based] → Python MAP[a-1:b] [0-based]
MAP = np.zeros(N)
MAP[0:1001]    = 42.0                          # Pre-CL (t=0-10 min)

# Transient response
MAP[1000:1401] = np.linspace(42, 84, 401)     # Phase 1: rise  step=0.105
MAP[1400:1801] = np.linspace(84, 73, 401)     # Phase 2: fall  step=0.0275

# First biased steady state — Group 2 attempt (will fail at 8.70 min)
MAP[1800:2401] = 73.0                          # Phase 3: t=18-24 min (6 min)

# Hemorrhage 1 (t=24-26 min)
MAP[2400:2601] = np.linspace(73, 57, 201)     # Hem1  step=-0.08

# Recovery to steady state
MAP[2600:3201] = np.linspace(57, 73, 601)     # Phase 4: step=0.02667

# Second biased steady state — part of Group 3 (will qualify)
MAP[3200:4201] = 73.0                          # Phase 5: t=32-42 min (10 min)

# Hemorrhage 2 (t=42-45 min)
MAP[4200:4501] = np.linspace(73, 57, 301)     # Hem2  step=-0.05333

# Second recovery
MAP[4500:5101] = np.linspace(57, 73, 601)     # Phase 6: step=0.02667

# Sinusoidal oscillation — non-zero Wobble
# Period=8 min, amplitude=4 mmHg, 2 complete cycles (t=51-67 min)
# MAP ∈ [69,77]: entirely within settling band [63,77]
# MAP=77 at peaks: PF=7/70=0.1 ≤ tol → in settling band but NOT in per10
t7 = Times[5100:6701]
MAP[5100:6701] = 73.0 + 4.0 * np.sin(2*np.pi/8 * (t7 - 51))

# Long post-oscillation steady state
MAP[6700:8501] = 73.0                          # Phase 8:   t=67-85 min
MAP[8500:8801] = 73.0                          # Phase 8.5: t=85-88 min

# Hemorrhage 3 (t=88-92 min) — near end → SSE Case 2
# Post-Hem3 = 8 min < 10 min → Case 2 triggered
MAP[8800:9201] = np.linspace(73, 55, 401)     # Hem3  step=-0.045

# Post-Hem3 recovery
MAP[9200:N]    = np.linspace(55, 70, N-9200)  # Phase 9: step=0.01875

# Hemorrhage signal
Hemorrhage = np.zeros(N)
Hemorrhage[2400:2601] = -0.1    # Hem1 (t=24-26 min)
Hemorrhage[4200:4501] = -0.1    # Hem2 (t=42-45 min)
Hemorrhage[8800:9201] = -0.1    # Hem3 (t=88-92 min)

# -----------------------------------------------------------------------
# 3. INDEPENDENT EXPECTED VALUE DERIVATIONS
# -----------------------------------------------------------------------
s       = ind_CL_start - 1    # 0-based: 1000
e       = ind_CL_end          # Python exclusive end: 10001
n_eng   = e - s               # 9001
PE_ind  = (MAP[s:e] - target) / target * 100.0
MAP_eng = MAP[s:e]

# (i) MDPE
# MAP>=73 samples: ~5006/9001 = 55.6% > 50% → median = 4.286%
MDPE_exp  = float(np.median(PE_ind))

# (ii) MDAPE
MDAPE_exp = float(np.median(np.abs(PE_ind)))

# (iii) Wobble
# Samples at exactly MAP=73: ~3410/9001 = 37.9% < 50% → Wobble > 0
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
# MAP_focus[0]=42, T10=44.8, T90=67.2, Phase 1 step=0.105
# idx_T10 (0-based local 27): MAP_focus[27]=42+27*0.105=44.835 → round=44.8 ✓
# idx_T90 (0-based local 240): MAP_focus[240]=42+240*0.105=67.2 → round=67.2 ✓
# risetime = (240-27)*0.01 = 2.13 min
MAP_focus = MAP[s:e]
T10_lvl   = round(MAP_focus[0] + 0.1*(target - MAP_focus[0]), 1)   # 44.8
T90_lvl   = round(MAP_focus[0] + 0.9*(target - MAP_focus[0]), 1)   # 67.2
r_map     = np.round(MAP_focus, 1)
idx_T10_exp = int(np.where(r_map >= T10_lvl)[0][0])   # local 0-based
idx_T90_exp = int(np.where(r_map >= T90_lvl)[0][0])   # local 0-based
risetime_exp = float((idx_T90_exp - idx_T10_exp) * Dt)  # 2.13 min

# (vii) Overshoot  max=84 → (84-70)/70*100 = 20%
overshoot_exp = float((np.max(MAP_focus) - target) / target * 100)

# (viii) Settling time — FOUR GROUPS
#
# Group 1 (Phase 1 upswing, local 0-based 200-333):
#   MAP[200]=42+200*0.105=63.0   PF=7/70=0.1 ≤ 0.1 ✓  (enters band)
#   MAP[334]=42+334*0.105=77.07  PF>0.1 ✗   MAP[333]=76.965 ✓ (exits)
#   Duration=(333-200)*0.01=1.34 min < 10 → FAIL
#
# Group 2 (Phase 2 downswing + Phase 3 + Hem1, local 655-1525):
#   Phase 2 step=-0.0275: MAP[655]=84-255*0.0275=76.9875 ≤ 77 ✓  (enters)
#   Hem1 exits at local 1525: MAP[1525]=73-125*0.08=63.0 ✓
#                              MAP[1526]=62.92 ✗
#   Duration=(1525-655)*0.01=8.70 min < 10 → FAIL (only 1.30 min short!)
#
# Group 3 (Phase 4+5+Hem2, local 1825-3387):
#   Phase 4 re-enters at local 1825: 57+225*(16/600)=63.0 ✓
#                                     MAP[1824]=62.973 ✗
#   Hem2 exits at local 3387: 73-187*(16/300)=63.027 ✓
#                              MAP[3388]=62.974 ✗
#   Duration=(3387-1825)*0.01=15.62 min ≥ 10 → PASS
#   settle = t_local[1825] = Times[s+1825] = Times[1000+1825] = Times[2825]
#   Times[2825] = 2825*0.01 = 28.25 min
#
# Group 4 (Phase 6+7+8+8.5+Hem3, local 3725-8022):
#   Phase 6 re-enters at local 3725 (same step as Phase 4)
#   Hem3 exits at local 8022: 73-222*0.045=63.01 ✓  MAP[8023]=62.965 ✗
#   Duration=(8022-3725)*0.01=42.97 min ≥ 10 → PASS (subsequent)
#
# MATLAB: settle = Times(ind_CL_start + 1826 - 1) = Times(2826) = 28.25 min
# Python: Times[ind_CL_start + 1826 - 2] = Times[2825] = 2825*0.01 = 28.25 min
settle_exp = float(Times[ind_CL_start + 1826 - 2])   # Times[2825] = 28.25 min

# (ix) SSE — Case 2 with three hemorrhage episodes
# Hem3 ends at global idx 9201 (t=92). CL ends at 10001 (t=100).
# Post-Hem3 = 8 min < 10 → Case 2.
# While loop walks back through all 3 episodes (~6500 iterations).
# Converges to before Hem1 onset (global ~2400, t≈24 min).
# 2-min window: Phase 3 (MAP=73) → SSE = 3 mmHg.
# Result holds regardless of convergence since ALL steady states = 73.
SSE_exp = 3.00

# -----------------------------------------------------------------------
# 4. RUN evaluate_CL_metrics
# -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end,
                              MAP, Times, target, Hemorrhage)

# -----------------------------------------------------------------------
# 5. COMPARISON TABLE
# -----------------------------------------------------------------------
tol = 1e-4

print()
print('=' * 67)
print('  Validation Script 3  —  Maximum Challenge Signal')
print('  3 hemorrhages | 4 settling groups | Wobble≠0 | SSE Case 2')
print('=' * 67)
print(f"  {'Metric':<22}  {'Expected':>10}  {'Computed':>10}  {'Status'}")
print('-' * 67)

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

print('-' * 67)
if all_pass:
    print(f'  All metrics matched within tolerance ({tol:.0e}).')
else:
    print('  One or more metrics did NOT match.')
print('=' * 67)

# -----------------------------------------------------------------------
# 6. DIAGNOSTICS
# -----------------------------------------------------------------------
print()
print('--- Settling group summary (threshold = 10.00 min) ---')
print(f"  {'Group':<6}  {'Start':>12}  {'End':>10}  {'Duration':>10}  Result")
groups = [
    (1, 200,  333,  1.34, 'FAIL (trivially short)'),
    (2, 655,  1525, 8.70, 'FAIL (only 1.30 min short of threshold)'),
    (3, 1825, 3387, 15.62,'PASS → settle = 28.25 min'),
    (4, 3725, 8022, 42.97,'PASS (subsequent, does not change result)'),
]
for g, start, end, dur, res in groups:
    print(f"  {g:<6}  {start:>12}  {end:>10}  {dur:>9.2f}m  {res}")

print()
print('--- Hemorrhage and SSE Case 2 summary ---')
print('  Hem1: t=24-26 min  (global idx 2400-2600)')
print('  Hem2: t=42-45 min  (global idx 4200-4500)')
print('  Hem3: t=88-92 min  (global idx 8800-9200)  ← last episode')
print('  Post-Hem3 to CL end: 8 min < 10 min → Case 2 triggered')
print('  While loop walks back through all 3 episodes (~6500 iterations)')
print('  Converges before Hem1 onset → window: Phase 3 (MAP=73)')
print(f'  SSE = {result["SSE"]:.2f} mmHg  (correct regardless of convergence)')

print()
print('--- PE indexing challenge ---')
PE_buggy   = (MAP[0:n_eng] - target) / target * 100.0
MDPE_buggy = float(np.median(PE_buggy))
print(f'  Correct MDPE (MAP[1000:10001]): {result["MDPE"]:7.3f} %')
print(f'  Buggy   MDPE (MAP[0:9001]):     {MDPE_buggy:7.3f} %')

print()
print('--- Wobble diagnostic ---')
pct_73 = float(np.sum(MAP[s:e] == 73.0)) / n_eng * 100
print(f'  Samples at MAP=73 exactly: {pct_73:.1f}%  (< 50%)')
print(f'  Phase 7 sinusoidal (period=8 min) spreads |PE-MDPE| away from 0')
print(f'  Wobble = {result["Wobble"]:.4f}%  (non-zero confirmed)')

# -----------------------------------------------------------------------
# 7. PLOT
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 6))
ax.set_xlim([0, 100])
ax.set_ylim([38, 92])

# Background shading
ax.axvspan(0,   10, alpha=0.10, color='gray')
ax.axvspan(10, 100, alpha=0.04, color='blue')
ax.axvspan(24,  26, alpha=0.25, color='red')
ax.axvspan(42,  45, alpha=0.25, color='red')
ax.axvspan(88,  92, alpha=0.30, color='red')

# MAP and reference lines
ax.plot(Times, MAP, 'b-', linewidth=1.5, label='MAP signal')
ax.axhline(y=target,     color='k', linestyle='--', linewidth=1.2,
           label='Target (70 mmHg)')
ax.axhline(y=target*1.1, color='r', linestyle=':',  linewidth=1.0)
ax.axhline(y=target*0.9, color='r', linestyle=':',  linewidth=1.0,
           label='±10% band')
ax.axvline(x=10, color='k', linestyle='-.', linewidth=1.0,
           label='CL start (t=10)')

# Rise time markers
ax.axvline(x=Times[s + idx_T10_exp], color='m', linestyle='--', linewidth=1.0)
ax.axvline(x=Times[s + idx_T90_exp], color='m', linestyle='--', linewidth=1.0)
ax.text(Times[s+idx_T10_exp]+0.4, T10_lvl+3, 'T₁₀', color='m', fontsize=8)
ax.text(Times[s+idx_T90_exp]+0.4, T90_lvl+3, 'T₉₀', color='m', fontsize=8)

# Settling time
ax.axvline(x=result['settlingtime'], color='g', linestyle='--',
           linewidth=1.2,
           label=f"Settle={result['settlingtime']:.2f} min (Group 3)")

# Group duration labels
ax.text(12, 88, 'G1:1.3m✗',  color='gray',              fontsize=7)
ax.text(18, 88, 'G2:8.7m✗',  color=(0.8, 0.4, 0.0),     fontsize=7)
ax.text(37, 88, 'G3:15.6m✓', color=(0.0, 0.6, 0.0),     fontsize=7)
ax.text(62, 88, 'G4:43.0m✓', color=(0.0, 0.4, 0.8),     fontsize=7)

# Metric annotations
ax.text(11, 83, f"Overshoot = {result['overshoot']:.1f}%",
        color='r', fontsize=8)
ax.text(11, 79, f"MDPE   = {result['MDPE']:.3f}%",
        color='b', fontsize=8)
ax.text(11, 75, f"Wobble = {result['Wobble']:.3f}%",
        color='b', fontsize=8)
ax.text(11, 71, f"SSE    = {result['SSE']:.2f} mmHg (Case 2)",
        color='b', fontsize=8)
ax.text(3, 58, 'Pre-CL', color='gray', fontsize=8, fontstyle='italic')

# Hemorrhage labels
for xt, label in [(24.1, 'Hem1'), (42.1, 'Hem2'), (88.1, 'Hem3')]:
    ax.text(xt, 40, label, color='darkred', fontsize=7)

# Legend
patch_pre = mpatches.Patch(color='gray', alpha=0.3, label='Pre-CL (excluded)')
patch_hem = mpatches.Patch(color='red',  alpha=0.3, label='Hemorrhages (×3)')
handles, leg_labels = ax.get_legend_handles_labels()
ax.legend(handles + [patch_pre, patch_hem],
          leg_labels + ['Pre-CL (excluded)', 'Hemorrhages (×3)'],
          loc='lower right', fontsize=7, ncol=2)

ax.set_xlabel('Time (min)')
ax.set_ylabel('MAP (mmHg)')
ax.set_title('Validation Script 3 — Maximum Challenge Signal\n'
             '3 Hemorrhages  |  4 Settling Groups  |  '
             'Sinusoidal Phase 7 (73±4 mmHg)  |  Wobble > 0')
ax.grid(True)
plt.tight_layout()
plt.show()
