%% =======================================================================
%  VALIDATION SCRIPT 1 FOR evaluate_CL_metrics
%
%  Generates a piecewise synthetic MAP signal with analytically known
%  metric values and verifies that evaluate_CL_metrics produces matching
%  results. Expected values are derived independently of the function
%  using direct formulas or closed-form reasoning.
%
%  Signal design (target = 70 mmHg, 100 Hz, 60 min):
%
%    Phase 1  t =  0 –  5 min   Linear rise  42 → 84 mmHg
%    Phase 2  t =  5 – 10 min   Linear fall  84 → 70 mmHg
%    Phase 3  t = 10 – 60 min   Constant     70 mmHg
%
%  ┌──────────────────────────────────────────────────────────────────┐
%  │  Analytically known values (derived below)                       │
%  │  Overshoot    = (84–70)/70 × 100        = 20.00 %               │
%  │  Rise time    = (idx_T90–idx_T10) × Dt  =  2.67 min             │
%  │  Settling     = Times(751)              =  7.50 min             │
%  │  SSE          = 70 – 70                 =  0.00 mmHg            │
%  │  MDPE/MDAPE/Wobble ≈ 0 % (Phase 3 dominates with 5001/6001      │
%  │                            samples at exactly the target)        │
%  │  Divergence   < 0 %/min  (|PE| decreases from ~40% to 0%)       │
%  └──────────────────────────────────────────────────────────────────┘
%
%  NOTE ON POT vs SETTLING BAND:
%    per10 uses strict inequality  |MAP–70| <  7 mmHg  (excludes MAP=77)
%    Settling uses ≤ inequality    |MAP–70| ≤  7 mmHg  (includes MAP=77)
%    This causes settling to begin one sample earlier (idx 751, t=7.50 min)
%    than the first sample counted by per10 (idx 752, t=7.51 min).
%
%  Reference: Chalumuri et al. Ann Biomed Eng 54, 857–868 (2026).
%% =======================================================================

clear; clc;

%% -----------------------------------------------------------------------
%  1. INPUTS
%% -----------------------------------------------------------------------
target       = 70;                     % mmHg
Dt           = 0.01;                   % min  (100 Hz)
Times        = (0 : Dt : 60);          % 1×6001 time vector (min)
N            = length(Times);          % 6001
ind_CL_start = 1;
ind_CL_end   = N;                      % full 60-min engagement
Hemorrhage   = zeros(1, N);            % no disturbance → SSE Case 1

%% -----------------------------------------------------------------------
%  2. SYNTHETIC SIGNAL CONSTRUCTION
%% -----------------------------------------------------------------------
%
%  Phase 1: linspace(42, 84, 501) → indices 1–501
%           MAP rises from 42 to 84 mmHg, step = 0.084 mmHg/sample
%  Phase 2: linspace(84, 70, 501) → indices 501–1001
%           MAP falls from 84 to 70 mmHg, step = 0.028 mmHg/sample
%           (MAP(501)=84 same value; overlap is harmless)
%  Phase 3: 70 mmHg constant     → indices 1001–6001
%           (MAP(1001)=70 same value; overlap is harmless)
%
MAP             = zeros(1, N);
MAP(1:501)      = linspace(42, 84, 501);
MAP(501:1001)   = linspace(84, 70, 501);
MAP(1001:N)     = 70;

%% -----------------------------------------------------------------------
%  3. INDEPENDENT EXPECTED VALUE DERIVATIONS
%% -----------------------------------------------------------------------

% ---- Performance Error over full engagement ----------------------------
PE_ind  = (MAP(ind_CL_start:ind_CL_end) - target) / target * 100;
n_eng   = ind_CL_end - ind_CL_start + 1;   % 6001 samples
MAP_eng = MAP(ind_CL_start:ind_CL_end);

% (i) MDPE [%]
% Phase 3 contributes 5001/6001 samples at PE = 0.
% The median of a distribution where the majority is at zero = 0%.
MDPE_exp = median(PE_ind);   % = 0.00 %

% (ii) MDAPE [%]
% Same reasoning as MDPE → 0.00 %
MDAPE_exp = median(abs(PE_ind));   % = 0.00 %

% (iii) Wobble [%]
% Wobble = median(|PE – MDPE|) = median(|PE – 0|) = median(|PE|) = 0.00 %
Wobble_exp = median(abs(PE_ind - MDPE_exp));   % = 0.00 %

% (iv) POT — percentage of time within X% of target [%]
% Uses strict inequality  |MAP – target| < X%×target
per10_exp = sum(abs(MAP_eng - target) < 0.1*target) * 100 / n_eng;
per20_exp = sum(abs(MAP_eng - target) < 0.2*target) * 100 / n_eng;
per30_exp = sum(abs(MAP_eng - target) < 0.3*target) * 100 / n_eng;
%
%  Manual derivation of per10 (for verification):
%    Phase 1 (idx 1–500):  step = 0.084 mmHg/sample
%      63 < MAP < 77  →  idx 252–417  →  166 samples
%    Phase 2 (idx 501–1000): step = 0.028 mmHg/sample
%      MAP < 77 at k > 751  →  idx 752–1000  →  249 samples
%    Phase 3 (idx 1001–6001):  MAP = 70  →  5001 samples
%    Total = 5416 / 6001 × 100 ≈ 90.25 %

% (v) Divergence [%/min]
% OLS slope of |PE| vs time.  |PE| starts at ~40% and decays to 0%;
% slope is expected to be negative (improving, stable control).
t_temp        = 0 : Dt : (length(PE_ind)-1)*Dt;
n_pe          = length(PE_ind);
divergence_exp = (sum(abs(PE_ind).*t_temp) - sum(abs(PE_ind))*mean(t_temp)) / ...
                 (sum(t_temp.^2) - sum(t_temp)^2/n_pe);
% Expected: negative value

% (vi) Rise time [min]
%  MAP_focus(1) = 42 mmHg   (< 0.9×70 = 63 → condition a met)
%  max(MAP)     = 84 mmHg   (> 70 → condition b met)
%
%  T10 = 42 + 0.1×(70–42) = 44.8 mmHg
%  T90 = 42 + 0.9×(70–42) = 67.2 mmHg
%
%  Phase 1 step = 0.084 mmHg/sample:
%    idx_T10: MAP(34) = 42+33×0.084 = 44.772 → round(44.772,1) = 44.8 ≥ 44.8 ✓
%             MAP(33) = 44.688         → round = 44.7 < 44.8                   ✗
%             → idx_T10 = 34
%    idx_T90: MAP(301) = 42+300×0.084 = 67.2  → round(67.2,1) = 67.2 ≥ 67.2  ✓
%             MAP(300) = 67.116         → round = 67.1 < 67.2                  ✗
%             → idx_T90 = 301
%
%  risetime = (301–34) × 0.01 = 2.67 min
MAP_focus   = MAP(ind_CL_start:ind_CL_end);
T10_lvl     = MAP_focus(1) + 0.1*(target - MAP_focus(1));   % 44.8
T90_lvl     = MAP_focus(1) + 0.9*(target - MAP_focus(1));   % 67.2
idx_T10_exp = min(find(round(MAP_focus,1) >= round(T10_lvl,1)));   % 34
idx_T90_exp = min(find(round(MAP_focus,1) >= round(T90_lvl,1)));   % 301
risetime_exp = (idx_T90_exp - idx_T10_exp) * Dt;                    % 2.67 min

% (vii) Overshoot [%]
%  max(MAP) = 84 mmHg → overshoot = (84–70)/70 × 100 = 20.00 %
overshoot_exp = (max(MAP_focus) - target) / target * 100;   % 20.00 %

% (viii) Settling time [min]
%  Settling uses |MAP–70|/70 ≤ 0.1 (note: ≤ not <, unlike per10)
%
%  Phase 1 group: MAP passes through [63,77] band during rise
%    MAP(251) = 63.0  → |63–70|/70 = 0.1 ≤ 0.1 ✓  (first in-band)
%    MAP(417) = 76.944 → in band ✓
%    MAP(418) = 77.028 → |77.028–70|/70 = 0.1004 > 0.1 ✗ (exits band)
%    Duration = Times(417)–Times(251) = 1.66 min < 10 min → does NOT qualify
%
%  Phase 2+3 group: MAP re-enters band as it falls through 77
%    MAP(751) = 84–250×0.028 = 77.0 → |77–70|/70 = 0.1 ≤ 0.1 ✓ (re-enters)
%    MAP(750) = 77.028 → 0.1004 > 0.1 ✗ (just outside)
%    From idx 751 onward MAP stays ≤ 77 and ≥ 70 → in band
%    Duration = Times(6001)–Times(751) = 52.50 min ≥ 10 min → QUALIFIES
%    → settling time = Times(751) = 7.50 min
settle_exp = Times(751);   % 7.50 min

% (ix) SSE [mmHg]
%  No hemorrhage → isempty(end_hem_ind) → Case 1:
%  Window = last 2 min = idx (6001–200+1):6001 = idx 5802:6001
%  MAP = 70 throughout Phase 3 → SSE = 70–70 = 0.00 mmHg
SSE_exp = mean(MAP(ind_CL_end - round(2/Dt) + 1 : ind_CL_end)) - target;   % 0.00

%% -----------------------------------------------------------------------
%  4. RUN evaluate_CL_metrics
%% -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end, ...
                              MAP, Times, target, Hemorrhage);

%% -----------------------------------------------------------------------
%  5. COMPARISON TABLE
%% -----------------------------------------------------------------------
tol = 1e-6;   % floating-point tolerance for PASS/FAIL

fprintf('\n');
fprintf('=================================================================\n');
fprintf('  Validation Results 1 for evaluate_CL_metrics\n');
fprintf('=================================================================\n');
fprintf('%-22s  %10s  %10s  %10s\n', 'Metric', 'Expected', 'Computed', 'Status');
fprintf('%s\n', repmat('-', 56, 1));

labels    = {'MDPE (%)',    'MDAPE (%)',   'Wobble (%)', ...
             'per10 (%)',   'per20 (%)',   'per30 (%)', ...
             'risetime (min)', 'overshoot (%)', 'settlingtime (min)', ...
             'divergence (%/min)', 'SSE (mmHg)'};

exp_vals  = [MDPE_exp,   MDAPE_exp,  Wobble_exp, ...
             per10_exp,  per20_exp,  per30_exp, ...
             risetime_exp, overshoot_exp, settle_exp, ...
             divergence_exp, SSE_exp];

comp_vals = [result.MDPE,    result.MDAPE,   result.Wobble, ...
             result.per10,   result.per20,   result.per30, ...
             result.risetime, result.overshoot, result.settlingtime, ...
             result.divergence, result.SSE];

all_pass = true;
for i = 1:length(labels)
    e = exp_vals(i);
    c = comp_vals(i);
    if abs(e - c) < tol
        status = 'PASS';
    else
        status = sprintf('FAIL (%.2e)', e - c);
        all_pass = false;
    end
    fprintf('%-22s  %10.4f  %10.4f  %s\n', labels{i}, e, c, status);
end

fprintf('%s\n', repmat('-', 56, 1));
if all_pass
    fprintf('  All metrics matched within tolerance (%.0e).\n', tol);
else
    fprintf('  One or more metrics did NOT match. Review output above.\n');
end
fprintf('=================================================================\n\n');

%% -----------------------------------------------------------------------
%  6. SIGNAL PLOT
%% -----------------------------------------------------------------------
figure; hold on; box on;
plot(Times, MAP, 'b-', 'LineWidth', 1.5);
yline(target,       'k--',  'LineWidth', 1.2);
yline(target*1.1,   'r:',   'LineWidth', 1.0);
yline(target*0.9,   'r:',   'LineWidth', 1.0);
xline(result.risetime*0 + Times(idx_T10_exp), 'm--', 'LineWidth', 1.0);
xline(result.risetime*0 + Times(idx_T90_exp), 'm--', 'LineWidth', 1.0);
xline(result.settlingtime, 'g--', 'LineWidth', 1.2);

text(Times(idx_T10_exp)+0.3, T10_lvl+2, 'T_{10}', 'Color','m','FontSize',9);
text(Times(idx_T90_exp)+0.3, T90_lvl+2, 'T_{90}', 'Color','m','FontSize',9);
text(result.settlingtime+0.3, 55, ...
     sprintf('Settle = %.2f min', result.settlingtime), ...
     'Color','g','FontSize',9);
text(31, 83, sprintf('Overshoot = %.1f%%', result.overshoot), ...
     'Color','r','FontSize',9);

xlabel('Time (min)');
ylabel('MAP (mmHg)');
title('Synthetic Validation Signal 1 — evaluate\_CL\_metrics');
legend('MAP signal', 'Target (70 mmHg)', '±10% band', ...
       'Location', 'southeast');
ylim([35 95]); grid on;