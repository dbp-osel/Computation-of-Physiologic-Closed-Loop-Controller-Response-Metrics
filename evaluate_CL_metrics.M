function compile = evaluate_CL_metrics(ind_CL_start, ind_CL_end, ...
MAP, Times, target, Hemorrhage)

%{
========================================================================
evaluate_CL_metrics - Computes  PCLC performance metrics
from a controlled physiologic variable response.
=======================================================================
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
------------------------------------------------------------------------

INPUTS:
------------------------------------------------------------------------
    ind_CL_start  - (scalar integer) Index in MAP and Times vectors
                    at which controller engagement begins.
                    Example: ind_CL_start = find(Times >= t_engage, 1)

    ind_CL_end    - (scalar integer) Index in MAP and Times vectors
                    at which controller engagement ends.
 
    MAP           - (vector, double) Full time-series of the controlled
                    physiologic variable (e.g., MAP in mmHg). Should cover 
		    the full duration of the test, including both pre-engagement 
		    and controller engagement periods.
                    Units: mmHg (for MAP-based applications) or units
                    of the controlled variable for other applications.

    Times         - (vector, double) Full time vector in minutes,
                    same length as MAP.
                    Units: minutes
 

    target        - (scalar, double) Target value of the controlled
                    physiologic variable that the controller is
                    regulating toward.
                    Example: target = 70 (mmHg) for MAP-based
                    fluid resuscitation.
                    Units: same as MAP.


    Hemorrhage    - (vector, double) Full time-series of the
                    disturbance signal (e.g., hemorrhage rate in
                    l/min for fluid resuscitation applications),
                    same length as MAP.
                    Used ONLY for determining the appropriate steady
                    state window in the SSE calculation. It is NOT
                    used in the computation of any other metric.
                    Set to a zero vector (zeros(size(MAP))) if no
                    disturbance was applied during the experiment.
                    Units: units of the disturbance signal
                    (e.g., l/min for hemorrhage rate).

------------------------------------------------------------------------

OUTPUTS:
------------------------------------------------------------------------

    compile       - (struct) Structure containing the following fields,
                    one value per subject:
        CV          Controlled Variable

        MDPE         Median Performance Error (%)
                     Quantifies systematic controller bias.
                     Positive: CV above target. Negative: CV below.
        MDAPE        Median Absolute Performance Error (%)
                     Quantifies overall controller accuracy.
        Wobble       Wobble (%)
                     Quantifies oscillatory behavior around target,
                     independent of systematic bias.
        per10        Percentage of time within 10% of target (%)
                     POT(10): proportion of engagement time where
                     |CV - target| < 0.1 * target.
        per20        Percentage of time within 20% of target (%)
                     POT(20): proportion of engagement time where
                     |CV - target| < 0.2 * target.
        per30        Percentage of time within 30% of target (%)
                     POT(30): proportion of engagement time where
                     |CV - target| < 0.3 * target.
        risetime     Modified rise time (min)
                     Time for CV to rise from 10% to 90% of the way
                     between CV at engagement and target.
                     Returns NaN if:
                       - CV at engagement >= 90% of target, OR
                       - CV never exceeds target during engagement.
        overshoot    Percentage overshoot (%)
                     Maximum percentage by which CV exceeds target.
                     Returns NaN if CV never exceeds target.
        settlingtime Settling time (min)
                     Earliest time CV enters and remains within 10%
                     of target for >= 10 continuous minutes.
                     Returns NaN if condition is never met.
        divergence   Divergence (%/min)
                     Slope of linear regression of |PE| over time.
                     Positive: increasing deviation from target.
                     Negative: decreasing deviation (stable control).
        SSE          Steady state error (mmHg, or units of CV)
                     Mean CV minus target during steady state window.
                     Positive: CV above target at steady state.
                     Negative: CV below target at steady state.
                     Window selection depends on disturbance signal
                     (see SSE section below for details).


REFERENCE:
------------------------------------------------------------------------
[1] Chalumuri et al. Annals of Biomedical Engineering, 54, pp.
    857-868 (2026). https://doi.org/10.1007/s10439-025-03929-2
[2] Varvel et al. Journal of Pharmacokinetics and
    Biopharmaceutics, 20(1), 63-94 (1992).
    https://doi.org/10.1007/bf01143186

-----------------------------------------------------------------------
Author: Ramin Bighamian, PhD
ramin.bighamian@fda.hhs.gov

For questions, contact: RST_CDRH@fda.hhs.gov
========================================================================
%}
%%-----------------------------------------------------------------------
%  SAMPLING INTERVAL
%  The tool is implemented at a fixed sampling rate of 100 Hz.
%  Dt = 0.01 min corresponds to 100 samples per minute.
%  If your data are differently sampled, revise the sampling rate accordingly.
%% -----------------------------------------------------------------------
Dt = 1/100; % minutes (100 Hz sampling rate)
%% -----------------------------------------------------------------------
%  PERFORMANCE ERROR (PE)
%  PE(t) = (CV(t) - target) / target * 100
%  Computed at each time instance during the controller engagement period.
%  PE is expressed as a percentage (%).
%% -----------------------------------------------------------------------
PE = (MAP(ind_CL_start:ind_CL_end) - target) / target * 100;


%% -----------------------------------------------------------------------
%  (i) MEDIAN PERFORMANCE ERROR (MDPE) [%]
%  Median of PE over the controller engagement period.
%  Quantifies systematic controller bias.
%  Positive MDPE: controlled variable systematically above target.
%  Negative MDPE: controlled variable systematically below target.
%  Reference: Varvel et al. [2], Chalumuri et al. [1]
%% -----------------------------------------------------------------------
MDPE = median(PE);
%% -----------------------------------------------------------------------
%  (ii) MEDIAN ABSOLUTE PERFORMANCE ERROR (MDAPE) [%]
%  Median of |PE| over the controller engagement period.
%  Measures overall controller accuracy without regard to the
%  direction of the error.
%  Reference: Varvel et al. [2], Chalumuri et al. [1]
%% -----------------------------------------------------------------------
MDAPE = median(abs(PE));
%% -----------------------------------------------------------------------
%  (iii) WOBBLE [%]
%  Median of |PE - MDPE| over the controller engagement period.
%  Quantifies oscillatory behavior of the controlled variable around
%  the target, independent of any systematic bias captured by MDPE.
%  Higher wobble indicates greater oscillation around the target.
%  Reference: Varvel et al. [2], Chalumuri et al. [1]

%% -----------------------------------------------------------------------

Wobble = median(abs(PE - MDPE));

%% -----------------------------------------------------------------------

%  (iv) PERCENTAGE OF TIME WITHIN X% OF TARGET [%]
%  POT(X): proportion of time during controller engagement that the
%  controlled variable remained within X% of the target value.
%  Three threshold levels are computed: 10%, 20%, and 30%.
%  For MAP target of 70 mmHg:
%    POT(10): [63, 77] mmHg
%    POT(20): [56, 84] mmHg
%    POT(30): [49, 91] mmHg
%  Note: uses the full engagement window ind_CL_start:ind_CL_end
%  (inclusive of both endpoints) for denominator computation.
%  Reference: Chalumuri et al. [1]

%% -----------------------------------------------------------------------
per10 = length(find(abs(MAP(ind_CL_start:ind_CL_end) - target) < ...
0.1 * target)) * 100 / (ind_CL_end - ind_CL_start + 1);

per20 = length(find(abs(MAP(ind_CL_start:ind_CL_end) - target) < ...
0.2 * target)) * 100 / (ind_CL_end - ind_CL_start + 1);

per30 = length(find(abs(MAP(ind_CL_start:ind_CL_end) - target) < ...
0.3 * target)) * 100 / (ind_CL_end - ind_CL_start + 1);

%% -----------------------------------------------------------------------
%  (v) DIVERGENCE [%/min]
%  Slope of the linear regression of |PE| over time during controller
%  engagement. Computed using an internal time vector constructed from
%  Dt, independent of the Times input, to align with PE vector indexing.
%  Positive divergence: increasing deviation from target over time
%                       (deteriorating control).
%  Negative divergence: decreasing deviation from target over time
%                       (improving, stable control).
%  Reference: Varvel et al. [2], Chalumuri et al. [1]

%% -----------------------------------------------------------------------
t_temp     = 0:Dt:(length(PE) - 1) * Dt;
divergence = (sum(abs(PE) .* t_temp) - ...
(sum(abs(PE))) * mean(t_temp)) / ...
(sum(t_temp .* t_temp) - ...
(sum(t_temp))^2 / length(t_temp));

%% -----------------------------------------------------------------------
%  (vi) MODIFIED RISE TIME [min]
%  Duration for the controlled variable to rise from 10% to 90% of
%  the way between its value at controller engagement (MAP_focus(1))
%  and the target.
%  Computed ONLY when BOTH of the following conditions are met:
%    (a) CV at engagement is below 90% of target (MAP_focus(1) < 0.9*target)
%    (b) CV exceeds the target at some point during engagement
%        (max(MAP_focus) > target)
%  If either condition is not met, risetime is returned as NaN.
%  Threshold crossings are evaluated on values rounded to 1 decimal
%  place to reduce sensitivity to noise near the threshold.
%  Scaling by Dt converts sample index difference to minutes.
%  Reference: Chalumuri et al. [1]
%% -----------------------------------------------------------------------
MAP_focus = MAP(ind_CL_start:ind_CL_end);

if MAP_focus(1) < 0.9 * target && max(MAP_focus) > target
    % Find first sample where CV reaches 10% of the way to target
    idx_T10 = min(find(round(MAP_focus, 1) >= ...
    round(MAP_focus(1) + 0.1 * (target - MAP_focus(1)), 1)));
    % Find first sample where CV reaches 90% of the way to target
    idx_T90 = min(find(round(MAP_focus, 1) >= ...
    round(MAP_focus(1) + 0.9 * (target - MAP_focus(1)), 1)));
    % Rise time in minutes
    risetime = (idx_T90 - idx_T10) * Dt;
else
    % Conditions for rise time computation not met
    risetime = nan;
end

%% -----------------------------------------------------------------------
%  (vii) PERCENTAGE OVERSHOOT [%]
%  Maximum percentage by which the controlled variable exceeds the
%  target during controller engagement.
%  Computed as: 100 * (max(CV) - target) / target
%  Returns NaN if the controlled variable never exceeds the target
%  during the engagement period (overshoot_aux <= 0).
%  Reference: Chalumuri et al. [1]
%% -----------------------------------------------------------------------
max_response  = max(MAP_focus);
overshoot_aux = 100 * (max_response - target) / target;

if overshoot_aux <= 0
    % CV never exceeded target: no overshoot occurred
    overshoot = nan;
else
    overshoot = overshoot_aux;
end
%% -----------------------------------------------------------------------
%  (viii) SETTLING TIME [min]
%  Earliest time at which the controlled variable enters and remains
%  within 10% of the target for a continuous duration of at least
%  10 minutes (settling_duration_threshold).
%  Algorithm:
%    1. Identify all time points where normalized absolute error
%       |CV - target| / target <= error_tol (10%)
%    2. Group these points into contiguous segments
%    3. Compute the duration of each segment
%    4. Return the start time of the first segment whose duration
%       meets or exceeds settling_duration_threshold (10 min)
%  Returns NaN if no qualifying segment exists.
%  Note: settling time is reported in the same units as Times (minutes).
%  Reference: Chalumuri et al. [1]

%% -----------------------------------------------------------------------
error_tol                   = 0.1;  % 10% normalized error tolerance
settling_duration_threshold = 10;   % minimum settling duration (minutes)
signal                      = MAP(ind_CL_start:ind_CL_end);
t                           = Times(ind_CL_start:ind_CL_end);
PF                          = abs(signal - target) / target;
settling_indices            = find(PF <= error_tol);
if ~isempty(settling_indices)
    % --- Identify start indices of each contiguous settling group ---
    % A new group starts when the previous sample was outside tolerance.
    % Always include the very first settling index as a group start,
    % then scan from i=2 to detect subsequent group starts.
    group_start_indices = [settling_indices(1)];
    for i = 2:length(settling_indices)
        current_index = settling_indices(i);
        past_index    = current_index - 1;
        if PF(past_index) > error_tol
            group_start_indices = [group_start_indices, current_index];
        end
    end

    % --- Identify end indices of each contiguous settling group ---
    % A group ends when the next sample goes outside tolerance.
    % Scan up to the second-to-last index, then append the very last
    % settling index as the end of the final group.
    group_end_indices = [];
    for i = 1:length(settling_indices) - 1
        current_index = settling_indices(i);
        next_index    = current_index + 1;
        if PF(next_index) > error_tol
            group_end_indices = [group_end_indices, current_index];
        end
    end
    % Always include the very last settling index as a group end
    group_end_indices = [group_end_indices, settling_indices(end)];

    % --- Compute duration of each settling group in minutes ---
    settling_durations = t(group_end_indices) - t(group_start_indices);

    % --- Retain only groups meeting the minimum duration threshold ---
    settling_groups = group_start_indices( ...
        settling_durations >= settling_duration_threshold);

    if ~isempty(settling_groups)
        % Return the start time of the first qualifying group
        start_time = t(settling_groups(1));
    else
        % No group met the minimum duration threshold
        start_time = nan;
    end
    settling = start_time;
else
    % CV never entered the tolerance band during engagement
    settling = nan;
end
%% -----------------------------------------------------------------------
%  (ix) STEADY STATE ERROR (SSE) [mmHg, or units of controlled variable]
%  Difference between mean controlled variable and target during the
%  steady state window.
%  Sign convention:
%    Positive SSE: CV above target at steady state
%    Negative SSE: CV below target at steady state
%
%  The metrics in this tool are designed for PCLC evaluation scenarios
%  subject to disturbances during controller engagement (e.g., hemorrhage
%  in fluid resuscitation). The disturbance signal (Hemorrhage) is used
%  SOLELY to determine the appropriate steady state window as follows:
%
%  WINDOW SELECTION LOGIC:
%    Case 1: If more than 10 minutes of data exist after the end of the
%            last disturbance, the system has had sufficient time to
%            reach steady state. The last 2 minutes of the engagement
%            period are used as the steady state window.
%    Case 2: If less than 10 minutes of data exist after the end of the
%            last disturbance, disturbances may have prevented the system
%            from reaching steady state by the end of the experiment.
%            The 2 minutes immediately before the last disturbance are
%            used as a proxy for steady state.
%
%  If no disturbance was applied, supply a zero vector for Hemorrhage.
%  In this case, end_hem_ind will be empty, and Case 1 will 
%  apply, using the last 2 minutes of the engagement period.
%
%  NOTE: A 20-point moving average smooth is applied to the Hemorrhage
%  signal internally ONLY for the purpose of disturbance end detection.
%  This does not affect the MAP data or any other metric computation.
%  Reference: Chalumuri et al. [1]
%% -----------------------------------------------------------------------
window            = 2;              % steady state window duration (min)
Hemorrhage_smooth = smooth(Hemorrhage, 20); % smooth disturbance signal
% for end-of-disturbance
% detection only

% Find the index of the last sample where disturbance is active
% (rounded smoothed Hemorrhage < 0 indicates active hemorrhage)
end_hem_ind  = max(find(round(Hemorrhage_smooth, 2) < 0));
last_hem_ind = length(Hemorrhage_smooth);


if isempty(end_hem_ind) || (last_hem_ind - end_hem_ind) * Dt > 10
    % Case 1: > 10 min after last disturbance
    % Use last 2 minutes of engagement period as steady state window
    sse = mean(MAP(ind_CL_end - round(window/Dt) + 1:ind_CL_end)) - target;
else
    % Case 2: < 10 min after last disturbance
    % Use 2 minutes before the last disturbance as steady state window
    % Walk back through the disturbance signal to find the onset of
    % the last disturbance episode

    % Initialize to the onset search start point before the loop
    cut_last_hem = end_hem_ind - 250;
    p = 0;
    while round(Hemorrhage_smooth( ...
        max(find(Hemorrhage_smooth(1:end_hem_ind - p - 250)))), ...
        2) ~= 0
        cut_last_hem = end_hem_ind - p - 250;
        p = p + 1;
    end
    start_idx = max(1, cut_last_hem - round(window/Dt) +1);  % guard against index < 1
    sse = mean(MAP(start_idx:cut_last_hem)) - target;
end



%% -----------------------------------------------------------------------

%  COMPILE OUTPUT STRUCT
%  All computed metrics are returned as fields of a single struct.
%  NaN is returned for metrics that could not be computed for this
%  subject (risetime, overshoot, settlingtime).
%  When reporting aggregate results across a cohort, NaN values should
%  be excluded from mean and standard deviation computations, and the
%  number of subjects for whom each metric was computable should be
%  reported.

%% -----------------------------------------------------------------------

compile = struct('MDPE',         MDPE,       ...
'MDAPE',        MDAPE,      ...
'Wobble',       Wobble,     ...
'per10',        per10,      ...
'per20',        per20,      ...
'per30',        per30,      ...
'risetime',     risetime,   ...
'overshoot',    overshoot,  ...
'settlingtime', settling,   ...
'divergence',   divergence, ...
'SSE',          sse);

end