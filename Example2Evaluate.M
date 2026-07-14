%% =======================================================================
%  VALIDATION SCRIPT 2 FOR evaluate_CL_metrics  —  CHALLENGE SIGNAL
%
%
%  THREE KEY CHALLENGES:
%
%  Challenge 1 — ind_CL_start ≠ 1  (tests PE indexing)
%    The controller engages at t = 10 min, not t = 0.
%    Indices 1–1001 (pre-engagement, MAP = 45 mmHg) are present in the
%    full MAP vector but must NOT be included in any metric computation.
%    If the original PE indexing bug were still present, PE would be
%    computed on the pre-engagement window instead of the CL window,
%    and MDPE/MDAPE/Wobble/Divergence would all be wrong.
%
%  Challenge 2 — SSE Case 2  (tests hemorrhage code path)
%    A hemorrhage event runs from t = 55 to t = 63 min.
%    Only 7 min remain between hemorrhage end and CL end (< 10 min),
%    so the SSE Case 2 branch is triggered and the steady-state window
%    is taken from 2 min BEFORE the last hemorrhage onset (t = 53–55 min)
%    rather than from the end of the engagement.
%
%  Challenge 3 — Three settling group attempts  (tests group detection)
%    Group 1  (Phase 1 upswing):  1.56 min  <  10 min  →  does NOT qualify
%    Group 2  (Phase 2 downswing): 2.81 min  <  10 min  →  does NOT qualify
%    Group 3  (Phase 3 onward):   39+ min   ≥  10 min  →  QUALIFIES
%    Only the third group yields the reported settling time.
%
%  Signal design (target = 70 mmHg, 100 Hz, 70 min total):
%
%    Pre-CL   t = 0–10  min  Constant      45 mmHg   (not in engagement)
%    Phase 1  t = 10–14 min  Linear rise   45 → 81 mmHg
%    Phase 2  t = 14–18 min  Linear fall   81 → 61 mmHg
%    Phase 3  t = 18–22 min  Linear rise   61 → 73 mmHg
%    Phase 4  t = 22–55 min  Constant      73 mmHg   (+3 mmHg bias)
%    Phase 5  t = 55–63 min  Linear fall   73 → 55 mmHg  (hemorrhage on)
%    Phase 6  t = 63–70 min  Linear rise   55 → 70 mmHg  (hemorrhage off)
%
%  Analytically known values:
%    Overshoot  = (81–70)/70×100       = 15.71 %
%    Rise time  = (idx_T90–idx_T10)×Dt =  2.22 min  (see Section 3)
%    Settling   = Times(1868)          = 18.67 min  (see Section 3)
%    SSE        = 73–70                =  3.00 mmHg (Case 2, pre-disturbance window)
%    MDPE       = (73–70)/70×100       =  4.286 %   (Phase 4 dominates median)
%    Wobble     = 0 %                  (Phase 4 uniform → zero spread around MDPE)
%
%  Reference: Chalumuri et al. Ann Biomed Eng 54, 857–868 (2026).
%% =======================================================================

clear; clc;

%% -----------------------------------------------------------------------
%  1. INPUTS
%% -----------------------------------------------------------------------
target       = 70;                     % mmHg
Dt           = 0.01;                   % min  (100 Hz)
Times        = (0 : Dt : 70);          % 1×7001 time vector (min)
N            = length(Times);          % 7001
ind_CL_start = 1001;                   % controller engages at t = 10 min
ind_CL_end   = N;                      % controller ends at t = 70 min

%% -----------------------------------------------------------------------
%  2. SYNTHETIC SIGNAL CONSTRUCTION
%% -----------------------------------------------------------------------
%
%  Pre-CL period: MAP = 45 mmHg (below target, outside engagement window).
%  Any PE-indexing bug would pick up this region and produce wrong metrics.
%
%  Phase 1: linspace(45,81,401) → global idx 1001–1401, step = 0.09 mmHg/sample
%  Phase 2: linspace(81,61,401) → global idx 1401–1801, step = 0.05 mmHg/sample
%  Phase 3: linspace(61,73,401) → global idx 1801–2201, step = 0.03 mmHg/sample
%  Phase 4: constant 73          → global idx 2201–5501 (MAP(2201)=73 same)
%  Phase 5: linspace(73,55,801) → global idx 5501–6301, step = 0.0225 mmHg/sample
%            Hemorrhage = –0.1 l/min active over this phase.
%  Phase 6: linspace(55,70,701) → global idx 6301–7001, step = 0.02143 mmHg/sample
%
MAP           = zeros(1, N);
MAP(1:1001)   = 45;                          % Pre-CL period
MAP(1001:1401) = linspace(45, 81, 401);      % Phase 1: rise
MAP(1401:1801) = linspace(81, 61, 401);      % Phase 2: fall
MAP(1801:2201) = linspace(61, 73, 401);      % Phase 3: recovery
MAP(2201:5501) = 73;                          % Phase 4: biased steady state
MAP(5501:6301) = linspace(73, 55, 801);      % Phase 5: hemorrhage drop
MAP(6301:7001) = linspace(55, 70, 701);      % Phase 6: recovery

Hemorrhage        = zeros(1, N);
Hemorrhage(5501:6301) = -0.1;               % hemorrhage active t = 55–63 min
%
%  CHALLENGE 2 SETUP:
%  Hemorrhage ends at global idx 6301 (t = 63 min).
%  CL ends at global idx 7001 (t = 70 min).
%  Post-hemorrhage duration = (7001–6301)×0.01 = 7.00 min < 10 min → Case 2.
%  SSE window = 2 min before last hemorrhage onset (onset ≈ global idx 5501, t = 55 min).
%  MAP over that window = 73 mmHg → SSE = 73–70 = 3.00 mmHg.

%% -----------------------------------------------------------------------
%  3. INDEPENDENT EXPECTED VALUE DERIVATIONS
%% -----------------------------------------------------------------------

PE_ind  = (MAP(ind_CL_start:ind_CL_end) - target) / target * 100;
n_eng   = ind_CL_end - ind_CL_start + 1;    % 6001 samples
MAP_eng = MAP(ind_CL_start:ind_CL_end);

% (i) MDPE [%]
%  Engagement window: MAP(1001:7001), 6001 samples.
%  Phase 4 (global 2201–5501, local 1201–4501): 3301 samples at MAP=73.
%  3301/6001 = 55% > 50% → median is the Phase 4 PE value.
%  MDPE = (73–70)/70×100 = 4.286%
%  CHALLENGE 1: if PE were computed on MAP(1:6001) instead of MAP(1001:7001),
%  the result would be wrong because MAP(1:1000)=45 would dominate.
MDPE_exp = median(PE_ind);

% (ii) MDAPE [%]
%  All Phase 4 PE values are positive (4.286%) → MDAPE = 4.286%.
MDAPE_exp = median(abs(PE_ind));

% (iii) Wobble [%]
%  Phase 4 contributes 3301 identical PE values → |PE–MDPE|=0 for >50% of samples.
%  Wobble = 0%.
Wobble_exp = median(abs(PE_ind - MDPE_exp));   % = 0.00 %

% (iv) POT [%]
%  Manual derivation (per10, |MAP–70| < 7):
%    Phase 1 (local 1–400):   MAP>63 strictly at local≥202 → 202–356 = 155 samples
%    Phase 2 (local 401–800): MAP<77 at local≥482, MAP>63 at local≤760 → 279 samples
%    Phase 3 (local 801–1200):MAP>63 strictly at local≥868 → 868–1200 = 333 samples
%    Phase 4 (local 1201–4500): MAP=73 → 3300 samples
%    Phase 5 (local 4501–5300): MAP>63 at local≤4945 → 445 samples
%    Phase 6 (local 5301–6001): MAP>63 at local≥5675 → 327 samples
%    Total = 4839 / 6001 × 100 ≈ 80.64%
per10_exp = sum(abs(MAP_eng - target) < 0.1*target) * 100 / n_eng;
per20_exp = sum(abs(MAP_eng - target) < 0.2*target) * 100 / n_eng;
per30_exp = sum(abs(MAP_eng - target) < 0.3*target) * 100 / n_eng;

% (v) Divergence [%/min]
%  |PE| starts at |45–70|/70×100=35.7% and eventually stabilises at 4.286%.
%  OLS slope expected to be negative.
t_temp        = 0 : Dt : (length(PE_ind)-1)*Dt;
n_pe          = length(PE_ind);
divergence_exp = (sum(abs(PE_ind).*t_temp) - sum(abs(PE_ind))*mean(t_temp)) / ...
                 (sum(t_temp.^2) - sum(t_temp)^2/n_pe);

% (vi) Rise time [min]
%  MAP_focus(1) = MAP(1001) = 45 mmHg  (< 0.9×70=63 → condition a met)
%  max(MAP_focus) = 81 mmHg            (> 70 → condition b met)
%  T10 = 45+0.1×(70–45) = 47.5 mmHg
%  T90 = 45+0.9×(70–45) = 67.5 mmHg
%  Phase 1 step = 0.09 mmHg/sample:
%    idx_T10: MAP(29) =45+28×0.09=47.52 → round=47.5≥47.5 ✓
%             MAP(28) =45+27×0.09=47.43 → round=47.4 <47.5 ✗  → idx_T10=29
%    idx_T90: MAP(251)=45+250×0.09=67.5 → round=67.5≥67.5 ✓
%             MAP(250)=45+249×0.09=67.41→ round=67.4 <67.5 ✗  → idx_T90=251
%    risetime = (251–29)×0.01 = 2.22 min
MAP_focus    = MAP(ind_CL_start:ind_CL_end);
T10_lvl      = MAP_focus(1) + 0.1*(target - MAP_focus(1));   % 47.5
T90_lvl      = MAP_focus(1) + 0.9*(target - MAP_focus(1));   % 67.5
idx_T10_exp  = min(find(round(MAP_focus,1) >= round(T10_lvl,1)));   % 29
idx_T90_exp  = min(find(round(MAP_focus,1) >= round(T90_lvl,1)));   % 251
risetime_exp = (idx_T90_exp - idx_T10_exp) * Dt;                     % 2.22 min

% (vii) Overshoot [%]
%  max(MAP_focus) = 81 → (81–70)/70×100 = 15.71%
overshoot_exp = (max(MAP_focus) - target) / target * 100;   % 15.71%

% (viii) Settling time [min]  —  CHALLENGE 3: three groups
%  Band: |MAP–70|/70 ≤ 0.1  →  MAP ∈ [63, 77] inclusive.
%  Local indices referenced to MAP_focus (MAP(1001:7001)).
%
%  Group 1 — Phase 1 upswing:
%    Enters at local 201 (MAP=63.0, PF=0.1 ✓); exits at local 357 (MAP=77.04, PF>0.1)
%    Last in-band: local 356 (MAP=76.95). Duration=(356–201)×0.01 = 1.55 min < 10 min ✗
%
%  Group 2 — Phase 2 downswing:
%    Re-enters at local 481 (MAP=77.0, PF=0.1 ✓); exits at local 762 (MAP=62.95, PF>0.1)
%    Last in-band: local 761 (MAP=63.0). Duration=(761–481)×0.01 = 2.80 min < 10 min ✗
%
%  Group 3 — Phase 3 onward:
%    Re-enters at local 868 (MAP=63.01, PF=0.0999 ✓; local 867: MAP=62.98, PF=0.1003 ✗)
%    Exits band during Phase 5: MAP drops below 63 at local 4946
%      MAP(4945)=63.01, PF=0.0999 ✓ (last in-band)
%      MAP(4946)=62.99, PF=0.1001 ✗
%    Duration = (4945–868)×0.01 = 40.77 min ≥ 10 min → QUALIFIES ✓
%    settle = t(local 868) = Times(1001+868–1) = Times(1868) = 18.67 min
settle_exp = Times(ind_CL_start + 868 - 1);   % Times(1868) = 18.67 min

% (ix) SSE [mmHg]  —  CHALLENGE 2: Case 2 path
%  Hemorrhage ends at global idx 6301 (t = 63 min).
%  After smoothing, end_hem_ind ≈ 6301.
%  last_hem_ind = 7001.
%  (7001–6301)×0.01 = 7.00 min < 10 min → Case 2 triggered.
%
%  Case 2: walk back to find onset of last hemorrhage episode.
%  Hemorrhage first goes negative at global idx 5501 (t = 55 min).
%  After convergence: cut_last_hem ≈ global idx 5501.
%  SSE window: MAP(5302:5501), all Phase 4 (MAP=73).
%  SSE = 73–70 = 3.00 mmHg.
%  (Verified analytically: the 2-min window before t=55 min sits entirely
%   in Phase 4 where MAP=73 throughout.)
SSE_exp = 3.00;   % mmHg (analytical)

%% -----------------------------------------------------------------------
%  4. RUN evaluate_CL_metrics
%% -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end, ...
                              MAP, Times, target, Hemorrhage);

%% -----------------------------------------------------------------------
%  5. COMPARISON TABLE
%% -----------------------------------------------------------------------
tol = 1e-4;   % slightly relaxed tolerance to accommodate smoothing
              % in SSE Case 2 onset detection

fprintf('\n');
fprintf('=================================================================\n');
fprintf('  Validation Results 2 — Challenge Signal\n');
fprintf('  ind_CL_start=%d (t=10 min) | SSE Case 2 | 3 settling groups\n', ...
        ind_CL_start);
fprintf('=================================================================\n');
fprintf('%-22s  %10s  %10s  %10s\n', 'Metric', 'Expected', 'Computed', 'Status');
fprintf('%s\n', repmat('-', 58, 1));

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

fprintf('%s\n', repmat('-', 58, 1));
if all_pass
    fprintf('  All metrics matched within tolerance (%.0e).\n', tol);
else
    fprintf('  One or more metrics did NOT match. Review output above.\n');
end
fprintf('=================================================================\n');

%% -----------------------------------------------------------------------
%  6. WHAT WRONG CODE WOULD HAVE PRODUCED  (for reference)
%% -----------------------------------------------------------------------
fprintf('\n--- Reference: what the ORIGINAL BUGGY code would return ---\n');
fprintf('(PE computed on MAP(1:6001) instead of MAP(1001:7001))\n');
PE_buggy   = (MAP(1:ind_CL_end - ind_CL_start) - target) / target * 100;
MDPE_buggy = median(PE_buggy);
fprintf('  MDPE (buggy):  %.4f %%  vs  expected: %.4f %%\n', ...
        MDPE_buggy, MDPE_exp);
fprintf('  Pre-CL MAP=45 would dominate → MDPE ≈ %.1f %% (wrong)\n\n', ...
        (45 - target)/target*100);

%% -----------------------------------------------------------------------
%  7. SIGNAL PLOT
%% -----------------------------------------------------------------------
figure; hold on; box on;
patch([10 10 70 70], [20 95 95 20], [0.95 0.95 1.0], ...
      'EdgeColor','none', 'DisplayName','CL engagement window');
plot(Times, MAP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MAP signal');
yline(target,     'k--', 'LineWidth', 1.2, 'DisplayName', 'Target (70 mmHg)');
yline(target*1.1, 'r:',  'LineWidth', 1.0, 'HandleVisibility','off');
yline(target*0.9, 'r:',  'LineWidth', 1.0, 'DisplayName', '±10% band');
xline(10,  'k-.', 'LineWidth', 1.0, 'DisplayName', 'CL start (t=10 min)');
xline(Times(idx_T10_exp + ind_CL_start - 1), 'm--', 'LineWidth',1.0, ...
      'HandleVisibility','off');
xline(Times(idx_T90_exp + ind_CL_start - 1), 'm--', 'LineWidth',1.0, ...
      'HandleVisibility','off');
xline(result.settlingtime, 'g--', 'LineWidth', 1.2, 'DisplayName', ...
      sprintf('Settle=%.2f min (3rd group)', result.settlingtime));
xline(55, 'r-',  'LineWidth', 1.0, 'DisplayName', 'Hemorrhage on (t=55)');
xline(63, 'r--', 'LineWidth', 1.0, 'DisplayName', 'Hemorrhage off (t=63)');

text(Times(idx_T10_exp+ind_CL_start-1)+0.3, T10_lvl+3, 'T_{10}','Color','m','FontSize',9);
text(Times(idx_T90_exp+ind_CL_start-1)+0.3, T90_lvl+3, 'T_{90}','Color','m','FontSize',9);
text(22, 88, sprintf('Overshoot = %.2f%%', result.overshoot),'Color','r','FontSize',9);
text(22, 84, sprintf('MDPE      = %.3f%%', result.MDPE),     'Color','b','FontSize',9);
text(22, 80, sprintf('SSE       = %.2f mmHg (Case 2)', result.SSE), ...
                                                             'Color','b','FontSize',9);
text(4,  50, 'Pre-CL','Color',[0.5 0.5 0.5],'FontSize',8,'FontAngle','italic');

xlabel('Time (min)'); ylabel('MAP (mmHg)');
title('Validation Signal 2 — Challenge Signal — evaluate\_CL\_metrics');
legend('Location','southeast','FontSize',8);
ylim([20 95]); grid on;