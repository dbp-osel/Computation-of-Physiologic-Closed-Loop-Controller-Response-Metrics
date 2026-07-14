"""
========================================================================
evaluate_CL_metrics.py
Computes PCLC performance metrics from a controlled physiologic
variable response.  Exact Python translation of evaluate_CL_metrics.m.
========================================================================

OVERVIEW:
This function computes a  set of closed-loop controller
performance metrics from time-series data of a controlled physiologic
variable collected during physiologic closed-loop controller (PCLC)
engagement. The function is presented as a representative example for
MAP-based PCLC evaluation in the context of automated fluid
resuscitation systems, where mean arterial pressure (MAP) is the
controlled variable and fluid infusion rate is the manipulated
variable. However, it can be applied to other physiologic variables
and PCLC applications.

The metrics are designed for PCLC evaluation scenarios in which the
controller may be subject to disturbances during the engagement
period. In the representative use case of automated fluid
resuscitation, these disturbances take the form of hemorrhagic events
that perturb MAP during closed-loop controller engagement. Several
metrics, including settling time and steady state error, include
specific provisions to account for the presence of such disturbances.
Users applying this function to other PCLC applications should supply
the appropriate disturbance signal for their specific context, or a
zero vector if no disturbance is applied.

NOTE ON DATA PREPARATION:
The function accepts the data as
provided and does not apply any additional processing to the
controlled variable signal. Data acquisition and any pre-processing applied can
be considered by the user prior to calling the function.

Disclaimer:

This software and documentation (the "Software") were developed at the Food
 and Drug Administration (FDA) by employees of the Federal Government in the
 course of their official duties. Pursuant to Title 17, Section 105 of the 
United States Code, this work is not subject to copyright protection and is
 in the public domain. Permission is hereby granted, free of charge, to any
 person obtaining a copy of the Software, to deal in the Software without 
restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, or sell copies of the Software or 
derivatives, and to permit persons to whom the Software is furnished to do 
so. FDA assumes no responsibility whatsoever for use by other parties of 
the Software, its source code, documentation or compiled executables, and 
makes no guarantees, expressed or implied, about its quality, reliability, 
or any other characteristic. Further, use of this code in no way implies 
endorsement by the FDA or confers any advantage in regulatory decisions. 
Although this software can be redistributed and/or modified freely, we ask 
that any derivative works bear some notice that they are derived from it, 
and any modified versions bear some notice that they have been modified.


NOTE ON INDEXING:
MATLAB uses 1-based indexing; Python uses 0-based indexing.
ind_CL_start and ind_CL_end are accepted here as 1-based integers
(same convention as the MATLAB function) and converted internally.
The rule used throughout is:
    MATLAB  MAP(ind_CL_start : ind_CL_end)   [1-based, both ends inclusive]
    Python  MAP[ind_CL_start-1 : ind_CL_end] [0-based, right end exclusive]

NOTE ON smooth():
MATLAB's smooth(y, 20) uses a 20-point moving average (span is
increased to 21 internally because it is even).  This is replicated
here with a centered rolling mean (min_periods=1 at boundaries).

REQUIREMENTS:
    pip install numpy pandas scipy matplotlib

REFERENCE:
[1] Chalumuri et al. Annals of Biomedical Engineering, 54, pp.
    857-868 (2026). https://doi.org/10.1007/s10439-025-03929-2
[2] Varvel et al. Journal of Pharmacokinetics and
    Biopharmaceutics, 20(1), 63-94 (1992).
    https://doi.org/10.1007/bf01143186

Author: Ramin Bighamian, PhD  (original MATLAB)
Python translation: direct port preserving all equations and logic.
========================================================================
"""

import math
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Internal helper: replicate MATLAB smooth(y, span) 'moving' average
# ---------------------------------------------------------------------------
def _matlab_smooth(y, span):
    """
    Centered moving average that replicates MATLAB's smooth(y, span).

    MATLAB increases an even span by 1, then averages span points
    centred at each position, shrinking the window at the edges
    (min_periods=1).
    """
    if span % 2 == 0:
        span += 1                    # MATLAB: even span → span + 1
    return (pd.Series(y)
              .rolling(window=span, center=True, min_periods=1)
              .mean()
              .values)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
def evaluate_CL_metrics(ind_CL_start, ind_CL_end,
                         MAP, Times, target, Hemorrhage):
    """
    Computes standardised PCLC performance metrics.

    Parameters
    ----------
    ind_CL_start : int
        1-based index at which controller engagement begins.
        Python equivalent of MATLAB: find(Times >= t_engage, 1)
        In Python:  np.argmax(Times >= t_engage) + 1
    ind_CL_end : int
        1-based index at which controller engagement ends.
        Python equivalent of MATLAB: length(Times)
        In Python:  len(Times)
    MAP : array-like
        Full time-series of the controlled physiologic variable
        (e.g. MAP in mmHg).
    Times : array-like
        Full time vector in minutes, same length as MAP.
    target : float
        Target value of the controlled physiologic variable.
    Hemorrhage : array-like
        Full time-series of the disturbance signal.
        Pass np.zeros(len(MAP)) if no disturbance was applied.

    Returns
    -------
    dict with keys:
        MDPE, MDAPE, Wobble,
        per10, per20, per30,
        risetime, overshoot, settlingtime,
        divergence, SSE
    All values are float.  NaN is returned where a metric cannot
    be computed (risetime, overshoot, settlingtime).
    """

    # ------------------------------------------------------------------
    # Input normalisation
    # ------------------------------------------------------------------
    MAP        = np.asarray(MAP,        dtype=float).ravel()
    Times      = np.asarray(Times,      dtype=float).ravel()
    Hemorrhage = np.asarray(Hemorrhage, dtype=float).ravel()

    # ------------------------------------------------------------------
    # SAMPLING INTERVAL
    # 100 Hz → Dt = 0.01 min
    # ------------------------------------------------------------------
    Dt = 1.0 / 100          # minutes

    # ------------------------------------------------------------------
    # Index conversion: MATLAB 1-based → Python 0-based slice notation
    #
    #   MATLAB: MAP(ind_CL_start : ind_CL_end)        inclusive both ends
    #   Python: MAP[ind_CL_start-1 : ind_CL_end]      0-based, excl. right
    #
    # s = ind_CL_start - 1  (0-based start)
    # e = ind_CL_end        (Python exclusive end = 1-based MATLAB end)
    # ------------------------------------------------------------------
    s = ind_CL_start - 1    # 0-based start
    e = ind_CL_end          # Python exclusive end

    # Number of engagement samples (= MATLAB ind_CL_end - ind_CL_start + 1)
    n_eng = e - s

    # ------------------------------------------------------------------
    # PERFORMANCE ERROR (PE)
    # PE(t) = (CV(t) - target) / target * 100   [%]
    # Computed over the full engagement window ind_CL_start:ind_CL_end.
    # ------------------------------------------------------------------
    PE = (MAP[s:e] - target) / target * 100.0

    # ------------------------------------------------------------------
    # (i) MEDIAN PERFORMANCE ERROR (MDPE)  [%]
    # Median of PE.  Positive → CV above target; negative → below.
    # ------------------------------------------------------------------
    MDPE = float(np.median(PE))

    # ------------------------------------------------------------------
    # (ii) MEDIAN ABSOLUTE PERFORMANCE ERROR (MDAPE)  [%]
    # Median of |PE|.  Measures accuracy irrespective of direction.
    # ------------------------------------------------------------------
    MDAPE = float(np.median(np.abs(PE)))

    # ------------------------------------------------------------------
    # (iii) WOBBLE  [%]
    # Median of |PE - MDPE|.  Oscillatory behaviour independent of bias.
    # ------------------------------------------------------------------
    Wobble = float(np.median(np.abs(PE - MDPE)))

    # ------------------------------------------------------------------
    # (iv) PERCENTAGE OF TIME WITHIN X% OF TARGET  [%]
    # Strict inequality |CV - target| < X * target, same as MATLAB.
    # Denominator = ind_CL_end - ind_CL_start + 1 (inclusive).
    # ------------------------------------------------------------------
    MAP_eng = MAP[s:e]

    per10 = float(np.sum(np.abs(MAP_eng - target) < 0.1 * target)
                  * 100.0 / n_eng)
    per20 = float(np.sum(np.abs(MAP_eng - target) < 0.2 * target)
                  * 100.0 / n_eng)
    per30 = float(np.sum(np.abs(MAP_eng - target) < 0.3 * target)
                  * 100.0 / n_eng)

    # ------------------------------------------------------------------
    # (v) DIVERGENCE  [%/min]
    # OLS slope of |PE| over time.
    # t_temp starts at 0 and increments by Dt, matching MATLAB's
    #   t_temp = 0:Dt:(length(PE)-1)*Dt
    # ------------------------------------------------------------------
    n_pe   = len(PE)
    t_temp = np.arange(n_pe) * Dt          # 0, Dt, 2Dt, ...
    absPE  = np.abs(PE)
    mean_t = float(np.mean(t_temp))

    divergence = float(
        (np.sum(absPE * t_temp) - np.sum(absPE) * mean_t)
        /
        (np.sum(t_temp ** 2) - np.sum(t_temp) ** 2 / n_pe)
    )

    # ------------------------------------------------------------------
    # (vi) MODIFIED RISE TIME  [min]
    # Duration for CV to rise from 10% to 90% of the way between
    # its value at engagement and the target.
    # Conditions (both must hold):
    #   (a) MAP_focus[0] < 0.9 * target
    #   (b) max(MAP_focus) > target
    # Threshold crossings evaluated on values rounded to 1 decimal place.
    # ------------------------------------------------------------------
    MAP_focus = MAP[s:e]           # engagement window (local 0-based)

    if MAP_focus[0] < 0.9 * target and float(np.max(MAP_focus)) > target:
        # Threshold levels (rounded to 1 d.p. as in MATLAB)
        T10_level = round(MAP_focus[0] + 0.1 * (target - MAP_focus[0]), 1)
        T90_level = round(MAP_focus[0] + 0.9 * (target - MAP_focus[0]), 1)

        rounded_map = np.round(MAP_focus, 1)

        arr_T10 = np.where(rounded_map >= T10_level)[0]
        arr_T90 = np.where(rounded_map >= T90_level)[0]

        if len(arr_T10) == 0 or len(arr_T90) == 0:
            risetime = float('nan')
        else:
            idx_T10  = int(arr_T10[0])     # first crossing (0-based local)
            idx_T90  = int(arr_T90[0])
            risetime = float((idx_T90 - idx_T10) * Dt)
    else:
        risetime = float('nan')

    # ------------------------------------------------------------------
    # (vii) PERCENTAGE OVERSHOOT  [%]
    # 100 * (max(CV) - target) / target.
    # NaN if CV never exceeds target.
    # ------------------------------------------------------------------
    max_response  = float(np.max(MAP_focus))
    overshoot_aux = 100.0 * (max_response - target) / target

    overshoot = float('nan') if overshoot_aux <= 0 else float(overshoot_aux)

    # ------------------------------------------------------------------
    # (viii) SETTLING TIME  [min]
    # Earliest time CV enters and stays within 10% of target for
    # >= 10 continuous minutes.
    #
    # All indices below are 0-based into the local engagement arrays
    # (signal, t_local, PF).  This is the direct Python equivalent of
    # MATLAB's 1-based indices into the same local arrays.
    # ------------------------------------------------------------------
    error_tol                   = 0.1    # 10% normalised error tolerance
    settling_duration_threshold = 10.0  # minimum duration [min]

    signal  = MAP[s:e]
    t_local = Times[s:e]
    PF      = np.abs(signal - target) / target

    settling_indices = np.where(PF <= error_tol)[0]    # 0-based

    if len(settling_indices) > 0:

        # --- Group start indices ---
        # Always include the first in-band sample as the first group start.
        # Then scan from position 1: a new group starts when the sample
        # immediately before the current one was outside tolerance.
        group_start_indices = [settling_indices[0]]
        for i in range(1, len(settling_indices)):
            current_index = settling_indices[i]
            past_index    = current_index - 1
            # past_index >= 0 is guaranteed because settling_indices[i] >= 1
            if PF[past_index] > error_tol:
                group_start_indices.append(current_index)
        group_start_indices = np.array(group_start_indices, dtype=int)

        # --- Group end indices ---
        # Scan up to second-to-last; a group ends when the next sample
        # is outside tolerance.  Always append the last in-band sample.
        group_end_indices = []
        for i in range(len(settling_indices) - 1):
            current_index = settling_indices[i]
            next_index    = current_index + 1
            if next_index < len(PF) and PF[next_index] > error_tol:
                group_end_indices.append(current_index)
        group_end_indices.append(int(settling_indices[-1]))
        group_end_indices = np.array(group_end_indices, dtype=int)

        # --- Durations [min] ---
        settling_durations = (t_local[group_end_indices]
                              - t_local[group_start_indices])

        # --- First qualifying group ---
        qualifying = group_start_indices[
            settling_durations >= settling_duration_threshold]

        settling = (float(t_local[qualifying[0]])
                    if len(qualifying) > 0 else float('nan'))
    else:
        settling = float('nan')

    # ------------------------------------------------------------------
    # (ix) STEADY STATE ERROR (SSE)  [units of CV]
    # Window selection:
    #   Case 1: > 10 min after last disturbance (or no disturbance)
    #           → last 2 min of engagement period.
    #   Case 2: < 10 min after last disturbance
    #           → 2 min immediately before last disturbance onset.
    #
    # 20-point moving-average smooth applied to Hemorrhage internally
    # only for disturbance-end detection (does not affect MAP or PE).
    #
    # INDEXING NOTE for SSE slices
    # ----------------------------
    # MATLAB (1-based):  MAP(a : b)   includes elements a through b.
    # Python  (0-based): MAP[a-1 : b] is the exact equivalent.
    #
    # Case 1:
    #   MATLAB:  MAP(ind_CL_end - w + 1 : ind_CL_end)   where w = 200
    #   Python:  MAP[ind_CL_end - w     : ind_CL_end]
    #
    # Case 2:
    #   cut_last_hem is stored as a 1-based integer (matching MATLAB).
    #   MATLAB:  MAP(start_idx : cut_last_hem)
    #            start_idx = max(1, cut_last_hem - w + 1)
    #   Python:  MAP[start_idx_0 : cut_last_hem]
    #            start_idx_0 = max(0, cut_last_hem - w)
    # ------------------------------------------------------------------
    window = 2.0                        # steady-state window duration [min]
    w      = round(window / Dt)         # = 200 samples

    # Smooth the disturbance signal (replicate MATLAB smooth(Hemorrhage,20))
    Hs = _matlab_smooth(Hemorrhage, 20)

    # Last sample where smoothed Hemorrhage rounds to < 0 (1-based style)
    active = np.where(np.round(Hs, 2) < 0)[0]   # 0-based indices

    last_hem_ind = len(Hs)              # = N (matches MATLAB length())

    if len(active) == 0:
        end_hem_ind = None              # no disturbance
    else:
        # Store as 1-based to keep all while-loop arithmetic identical to MATLAB
        end_hem_ind = int(active[-1]) + 1

    if end_hem_ind is None or (last_hem_ind - end_hem_ind) * Dt > 10:
        # ---- Case 1 ----
        # Python: MAP[ind_CL_end - w : ind_CL_end]
        #         = MATLAB MAP(ind_CL_end - w + 1 : ind_CL_end)  ✓ (200 samples)
        sse = float(np.mean(MAP[ind_CL_end - w : ind_CL_end]) - target)

    else:
        # ---- Case 2 ----
        # Walk back to find onset of the last disturbance episode.
        # cut_last_hem is maintained as a 1-based integer throughout,
        # matching the MATLAB variable of the same name exactly.
        cut_last_hem = end_hem_ind - 250   # 1-based initialisation
        p = 0

        while True:
            window_end = end_hem_ind - p - 250   # 1-based; = length of subarray

            # Exit when window is empty (matches MATLAB while([]) → false)
            if window_end <= 0:
                break

            # MATLAB: Hemorrhage_smooth(1 : window_end)
            # Python: Hs[0 : window_end]  (first window_end elements)
            subarray = Hs[:window_end]

            # MATLAB: max(find(subarray))  → index of last nonzero element
            nonzero_idx = np.where(subarray)[0]   # 0-based

            # Exit when no nonzero found (matches MATLAB max([]) → empty)
            if len(nonzero_idx) == 0:
                break

            last_val = subarray[nonzero_idx[-1]]

            # Exit when last nonzero rounds to 0
            if np.round(last_val, 2) == 0:
                break

            # Update (same assignments as MATLAB loop body)
            cut_last_hem = end_hem_ind - p - 250
            p += 1

        # MATLAB: start_idx = max(1, cut_last_hem - w + 1)   [1-based]
        # Python: start_idx_0 = max(0, cut_last_hem - w)      [0-based]
        start_idx_0 = max(0, cut_last_hem - w)

        # MATLAB: sse = mean(MAP(start_idx : cut_last_hem)) - target
        # Python: MAP[start_idx_0 : cut_last_hem]
        #   start_idx_0 is 0-based start;
        #   cut_last_hem is 1-based end = Python exclusive end  ✓
        sse = float(np.mean(MAP[start_idx_0 : cut_last_hem]) - target)

    # ------------------------------------------------------------------
    # COMPILE OUTPUT
    # ------------------------------------------------------------------
    return {
        'MDPE':         MDPE,
        'MDAPE':        MDAPE,
        'Wobble':       Wobble,
        'per10':        per10,
        'per20':        per20,
        'per30':        per30,
        'risetime':     risetime,
        'overshoot':    overshoot,
        'settlingtime': settling,
        'divergence':   divergence,
        'SSE':          sse,
    }
