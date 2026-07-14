%% =======================================================================
%  VALIDATION SCRIPT 3  —  MAXIMUM CHALLENGE SIGNAL
%
%  The most comprehensive validation signal, designed to stress-test
%  every metric simultaneously with features not combined in any single
%  previous script.
%
%  =======================================================================
%  WHAT IS NEW / DIFFERENT FROM EARLIER SCRIPTS
%  =======================================================================
%
%  Script 1–2: ind_CL_start=1, no hemorrhage, constant steady state.
%  Script 3:   ind_CL_start≠1, SSE Case 2, 3 settling groups.
%  Script 4:   ind_CL_start≠1, 2 hemorrhages, 3 settling groups, Wobble≠0.
%
%  This script adds:
%    • THREE hemorrhage events (not one or two)
%    • FOUR settling groups (not three) — two fail, two qualify
%    • Group 2 fails at 8.70 min (dramatically close to 10 min threshold)
%    • Group 3 qualifies at 15.62 min (clearly above threshold)
%    • Sinusoidal oscillation that drives Wobble non-zero
%    • All three hemorrhage episodes interact with the SSE Case 2 logic
%    • ind_CL_start = 1001 (pre-CL period at MAP = 42, same as CL start
%      value, so only the engagement window gives the correct MDPE)
%
%  =======================================================================
%  SIGNAL DESIGN  (target = 70 mmHg, 100 Hz, 100 min total)
%  =======================================================================
%
%    Pre-CL    t =  0– 10 min   Constant        42 mmHg  (not engaged)
%    Phase 1   t = 10– 14 min   Linear rise     42 → 84 mmHg
%    Phase 2   t = 14– 18 min   Linear fall     84 → 73 mmHg
%    Phase 3   t = 18– 24 min   Constant        73 mmHg   (Group 2 attempt)
%    Hem1      t = 24– 26 min   Linear fall     73 → 57 mmHg
%    Phase 4   t = 26– 32 min   Linear rise     57 → 73 mmHg
%    Phase 5   t = 32– 42 min   Constant        73 mmHg
%    Hem2      t = 42– 45 min   Linear fall     73 → 57 mmHg
%    Phase 6   t = 45– 51 min   Linear rise     57 → 73 mmHg
%    Phase 7   t = 51– 67 min   Sinusoidal      73 ± 4 mmHg  (Wobble≠0)
%    Phase 8   t = 67– 85 min   Constant        73 mmHg
%    Phase 8.5 t = 85– 88 min   Constant        73 mmHg   (SSE window)
%    Hem3      t = 88– 92 min   Linear fall     73 → 55 mmHg
%    Phase 9   t = 92–100 min   Linear rise     55 → 70 mmHg
%
%  =======================================================================
%  FOUR SETTLING GROUPS  (|MAP−70|/70 ≤ 0.1,  i.e.  MAP ∈ [63, 77])
%  =======================================================================
%
%    Group 1  Phase 1 upswing          local 201–334     1.34 min  ✗ FAIL
%    Group 2  Phase 2+3+Hem1           local 656–1526    8.70 min  ✗ FAIL
%             (barely fails — only 1.3 min short of threshold)
%    Group 3  Phase 4+5+Hem2           local 1826–3388  15.62 min  ✓ PASS
%             → settle = Times(2826) = 28.25 min
%    Group 4  Phase 6+7+8+8.5+Hem3    local 3726–8023  42.97 min  ✓ PASS
%             (subsequent; settle already found in Group 3)
%
%  =======================================================================
%  ANALYTICALLY KNOWN EXPECTED VALUES
%  =======================================================================
%
%    MDPE        = (73−70)/70×100          =  4.286 %
%                  Confirmed: 55.6% of CL samples have MAP ≥ 73
%    Wobble      > 0 %  (sinusoidal Phase 7: <42% of samples at MAP=73)
%    Overshoot   = (84−70)/70×100          = 20.000 %
%    Rise time   = (idx_T90−idx_T10)×Dt   =  2.13 min
%    Settle time = Times(2826)             = 28.25 min  (Group 3)
%    SSE         = 73−70                   =  3.000 mmHg  (Case 2)
%                  All pre-hemorrhage steady states at 73 → window always
%                  returns 73 regardless of which episode onset is found
%
%  Reference: Chalumuri et al. Ann Biomed Eng 54, 857–868 (2026).
%% =======================================================================

clear; clc;

%% -----------------------------------------------------------------------
%  1. INPUTS
%% -----------------------------------------------------------------------
target       = 70;
Dt           = 0.01;
Times        = (0 : Dt : 100);         % 1×10001
N            = length(Times);          % 10001
ind_CL_start = 1001;                   % t = 10 min
ind_CL_end   = N;                      % t = 100 min

%% -----------------------------------------------------------------------
%  2. SIGNAL CONSTRUCTION
%% -----------------------------------------------------------------------
MAP = zeros(1, N);
MAP(1:1001) = 42;                                    % Pre-CL (t=0–10 min)

% ---- Transient response ------------------------------------------------
MAP(1001:1401) = linspace(42, 84, 401);  % Phase 1: rise  step=0.105
MAP(1401:1801) = linspace(84, 73, 401);  % Phase 2: fall  step=0.0275

% ---- First biased steady state  (Group 2 attempt, will fail) -----------
% Duration of Group 2 from Phase 2 band-entry to Hem1 band-exit:
% t = 16.55 → 25.26 min = 8.70 min  <  10 min threshold → FAIL
MAP(1801:2401) = 73;                     % Phase 3: t=18–24 min (6 min)

% ---- Hemorrhage 1  (t = 24–26 min) ------------------------------------
MAP(2401:2601) = linspace(73, 57, 201);  % Hem1  step = –0.08

% ---- Recovery to biased steady state -----------------------------------
MAP(2601:3201) = linspace(57, 73, 601);  % Phase 4: t=26–32  step=0.02667

% ---- Second biased steady state  (part of Group 3, which qualifies) ----
MAP(3201:4201) = 73;                     % Phase 5: t=32–42 min (10 min)

% ---- Hemorrhage 2  (t = 42–45 min) -----------------------------------
MAP(4201:4501) = linspace(73, 57, 301);  % Hem2  step = –0.05333

% ---- Second recovery ---------------------------------------------------
MAP(4501:5101) = linspace(57, 73, 601);  % Phase 6: t=45–51  step=0.02667

% ---- Sinusoidal oscillation  (introduces non-zero Wobble) --------------
%  Period = 8 min = 800 samples, amplitude = 4 mmHg, 2 complete cycles
%  MAP ∈ [69, 77] — entirely within settling band [63, 77]
%  MAP = 77 at peaks (sin=1): PF = 7/70 = 0.1 ≤ tol → still in band
%  MAP = 77 excluded by per10 (strict < 7) but included by settling (≤ 7)
t7 = Times(5101:6701);
MAP(5101:6701) = 73 + 4 * sin(2*pi/8 * (t7 - 51));  % Phase 7: t=51–67

% ---- Long post-oscillation steady state --------------------------------
MAP(6701:8501) = 73;                     % Phase 8:   t=67–85 min (18 min)
MAP(8501:8801) = 73;                     % Phase 8.5: t=85–88 min (SSE buffer)

% ---- Hemorrhage 3  (t = 88–92 min) — near end → SSE Case 2 -----------
%  Hem3 ends at t=92. CL ends at t=100. Post-Hem3 = 8 min < 10 → Case 2.
%  With 3 episodes, the while loop converges to before Hem1 (t≈24).
%  2-min window before Hem1 onset = Phase 3 (MAP=73). SSE = 3 mmHg.
MAP(8801:9201) = linspace(73, 55, 401);  % Hem3  step = –0.045

% ---- Post-Hem3 recovery ------------------------------------------------
MAP(9201:N)    = linspace(55, 70, 801);  % Phase 9: t=92–100  step=0.01875

% ---- Hemorrhage signal -------------------------------------------------
Hemorrhage            = zeros(1, N);
Hemorrhage(2401:2601) = -0.1;            % Hem1 (t = 24–26 min)
Hemorrhage(4201:4501) = -0.1;            % Hem2 (t = 42–45 min)
Hemorrhage(8801:9201) = -0.1;            % Hem3 (t = 88–92 min)

%% -----------------------------------------------------------------------
%  3. INDEPENDENT EXPECTED VALUE DERIVATIONS
%% -----------------------------------------------------------------------

PE_ind  = (MAP(ind_CL_start:ind_CL_end) - target) / target * 100;
n_eng   = ind_CL_end - ind_CL_start + 1;    % 9001 local samples
MAP_eng = MAP(ind_CL_start:ind_CL_end);

% (i) MDPE [%]
%  Samples with MAP≥73 in CL engagement window (local 1–9001):
%    Phase 1 (local 297–400): 104 samples MAP∈(73,84)
%    Phase 2 non-overlap (local 401–800): 400 samples MAP∈(73,84)
%    Phase 3 non-overlap (local 801–1400): 600 @ MAP=73
%    Phase 5 non-overlap (local 2201–3200): 1000 @ MAP=73
%    Phase 7 positive half-cycles: ~800 samples MAP∈(73,77]
%    Phase 8 non-overlap (local 5701–7500): 1800 @ MAP=73
%    Phase 8.5 non-overlap (local 7501–7800): 300 @ MAP=73
%    Hem2/Hem3 first-sample endpoints: ~2 @ MAP=73
%  Total ≈ 5006 / 9001 = 55.6% > 50% → median = 4.286%
MDPE_exp = median(PE_ind);   % 4.286%

% (ii) MDAPE [%]
%  Phase 3/5/8/8.5 all at MAP=73 (PE=4.286%, |PE|=4.286%).
%  Phase 4/6 rise through (57→73): most samples below 73 → |PE|>4.286%.
%  Independent computation required.
MDAPE_exp = median(abs(PE_ind));

% (iii) Wobble [%]
%  Samples at exactly MAP=73: ~3717/9001 = 41.3% < 50%.
%  Phase 7 sinusoidal spreads |PE–MDPE| = |4sin(θ)|/70×100 ∈ [0, 5.71%].
%  Wobble > 0 guaranteed since <50% of samples contribute zero.
Wobble_exp = median(abs(PE_ind - MDPE_exp));

% (iv) POT [%]
%  per10 uses strict < 7 mmHg; per20/30 analogous.
%  Note: Phase 7 peaks at MAP=77 → |77−70|=7 NOT < 7 → excluded from per10.
per10_exp = sum(abs(MAP_eng - target) < 0.1*target) * 100 / n_eng;
per20_exp = sum(abs(MAP_eng - target) < 0.2*target) * 100 / n_eng;
per30_exp = sum(abs(MAP_eng - target) < 0.3*target) * 100 / n_eng;

% (v) Divergence [%/min]
%  |PE| starts at (42−70)/70×100 = 40% and declines to ~4.286%.
%  Expected negative slope (improving control).
t_temp        = 0 : Dt : (length(PE_ind)-1)*Dt;
n_pe          = length(PE_ind);
divergence_exp = (sum(abs(PE_ind).*t_temp) - sum(abs(PE_ind))*mean(t_temp)) / ...
                 (sum(t_temp.^2) - sum(t_temp)^2/n_pe);

% (vi) Rise time [min]
%  MAP_focus(1) = 42 mmHg  (< 0.9×70=63  ✓  condition a)
%  max(MAP_focus) = 84 mmHg  (> 70  ✓  condition b)
%
%  T10 = 42 + 0.1×(70−42) = 44.8 mmHg
%  T90 = 42 + 0.9×(70−42) = 67.2 mmHg
%
%  Phase 1 step = (84−42)/400 = 0.105 mmHg/sample:
%    idx_T10: MAP(28) = 42+27×0.105 = 44.835 → round(44.835,1) = 44.8 ≥ 44.8 ✓
%             MAP(27) = 44.73              → round = 44.7 < 44.8            ✗
%             → idx_T10 = 28
%    idx_T90: MAP(241)= 42+240×0.105 = 67.2  → round(67.2,1) = 67.2 ≥ 67.2 ✓
%             MAP(240)= 67.095             → round = 67.1 < 67.2            ✗
%             → idx_T90 = 241
%
%    risetime = (241−28)×0.01 = 213×0.01 = 2.13 min
MAP_focus    = MAP(ind_CL_start:ind_CL_end);
T10_lvl      = MAP_focus(1) + 0.1*(target - MAP_focus(1));   % 44.8
T90_lvl      = MAP_focus(1) + 0.9*(target - MAP_focus(1));   % 67.2
idx_T10_exp  = min(find(round(MAP_focus,1) >= round(T10_lvl,1)));   % 28
idx_T90_exp  = min(find(round(MAP_focus,1) >= round(T90_lvl,1)));   % 241
risetime_exp = (idx_T90_exp - idx_T10_exp) * Dt;                     % 2.13 min

% (vii) Overshoot [%]
%  max(MAP_focus) = 84  →  (84−70)/70×100 = 20.000%
overshoot_exp = (max(MAP_focus) - target) / target * 100;

% (viii) Settling time [min]  —  FOUR GROUPS
%
%  GROUP 1 — Phase 1 upswing  (trivially short):
%    Phase 1 step = 0.105. MAP∈[63,77]:
%    Enters: MAP(201) = 42+200×0.105 = 63.0 → PF=7/70=0.1 ≤ 0.1 ✓
%            MAP(200) = 62.895              → PF>0.1                   ✗
%    Exits:  MAP(335) = 42+334×0.105 = 77.07 → PF>0.1                 ✗
%            MAP(334) = 76.965              → PF≤0.1                   ✓ (last in)
%    Duration = (334−201)×0.01 = 1.34 min  < 10 min  → FAIL ✗
%
%  GROUP 2 — Phase 2 downswing + Phase 3 + into Hem1  (barely fails):
%    Phase 2 step = –(84–73)/400 = –0.0275. MAP≤77:
%    Enters: MAP(656) = 84−255×0.0275 = 76.9875 → PF=6.9875/70≤0.1  ✓
%            MAP(655) = 77.015               → PF>0.1                  ✗
%    Phase 3 constant 73 (all in band).
%    Hem1: linspace(73,57,201), step=–0.08. MAP<63:
%    MAP(1527) = 73−126×0.08 = 62.92 → PF>0.1                        ✗
%    MAP(1526) = 73−125×0.08 = 63.0  → PF=7/70=0.1 ≤ 0.1            ✓ (last in)
%    Duration = (1526−656)×0.01 = 8.70 min  < 10 min  → FAIL ✗
%    (Only 1.30 min short of the threshold — the most demanding near-miss)
%
%  GROUP 3 — Phase 4 recovery + Phase 5 + into Hem2  (qualifies):
%    Phase 4: linspace(57,73,601), step=16/600≈0.02667. MAP≥63:
%    Enters: MAP(1826) = 57+225×(16/600) = 63.0 → PF=0.1 ≤ 0.1      ✓
%            MAP(1825) = 62.973                  → PF>0.1              ✗
%    Phase 5 constant 73 (all in band).
%    Hem2: linspace(73,57,301), step=–16/300. MAP<63:
%    MAP(3389) = 73−188×(16/300) = 62.973 → PF>0.1                   ✗
%    MAP(3388) = 73−187×(16/300) = 63.027 → PF≤0.1                   ✓ (last in)
%    Duration = (3388−1826)×0.01 = 15.62 min  ≥ 10 min  → PASS ✓
%    settle = t(local 1826) = Times(1001+1826–1) = Times(2826) = 28.25 min
%
%  GROUP 4 — Phase 6+7+8+8.5+into Hem3  (long subsequent group):
%    Phase 6 re-enters at:
%    MAP(3726) = 57+225×(16/600) = 63.0 → ✓; MAP(3725)=62.973 → ✗
%    Hem3 exits at MAP(8024)=62.965 → ✗; MAP(8023)=63.01 → ✓
%    Duration = (8023−3726)×0.01 = 42.97 min  ≥ 10 min  → PASS ✓
%    (settle already found in Group 3; Group 4 does not change result)
settle_exp = Times(ind_CL_start + 1826 - 1);   % Times(2826) = 28.25 min

% (ix) SSE [mmHg]  —  THREE HEMORRHAGE EPISODES, Case 2
%
%  Hem3 ends at global 9201 (t=92). After smoothing, end_hem_ind≈9211.
%  last_hem_ind = 10001.
%  (10001–9211)×0.01 = 7.9 min  < 10 min  →  Case 2 triggered.
%
%  The while loop walks backward through Hem3, Hem2, and Hem1 sequentially
%  before the search window clears all hemorrhage (requires ~6500+ iterations).
%  The loop converges to a position before Hem1 onset (global 2401, t=24 min).
%  2-min window before Hem1 onset: global 2202–2401 → Phase 3 (MAP=73).
%  SSE = 73–70 = 3.00 mmHg.
%
%  This result holds regardless of which episode onset the loop finds
%  because ALL pre-hemorrhage steady-state phases have MAP=73 mmHg.
SSE_exp = 3.00;

%% -----------------------------------------------------------------------
%  4. RUN evaluate_CL_metrics
%% -----------------------------------------------------------------------
result = evaluate_CL_metrics(ind_CL_start, ind_CL_end, ...
                              MAP, Times, target, Hemorrhage);

%% -----------------------------------------------------------------------
%  5. COMPARISON TABLE
%% -----------------------------------------------------------------------
tol = 1e-4;

fprintf('\n');
fprintf('=================================================================\n');
fprintf('  Validation Script 3  —  Maximum Challenge Signal\n');
fprintf('  3 hemorrhages | 4 settling groups | Wobble≠0 | SSE Case 2\n');
fprintf('=================================================================\n');
fprintf('%-22s  %10s  %10s  %10s\n','Metric','Expected','Computed','Status');
fprintf('%s\n', repmat('-',58,1));

labels = {'MDPE (%)', 'MDAPE (%)', 'Wobble (%)', ...
          'per10 (%)', 'per20 (%)', 'per30 (%)', ...
          'risetime (min)', 'overshoot (%)', 'settlingtime (min)', ...
          'divergence (%/min)', 'SSE (mmHg)'};

exp_vals  = [MDPE_exp, MDAPE_exp, Wobble_exp, ...
             per10_exp, per20_exp, per30_exp, ...
             risetime_exp, overshoot_exp, settle_exp, ...
             divergence_exp, SSE_exp];

comp_vals = [result.MDPE, result.MDAPE, result.Wobble, ...
             result.per10, result.per20, result.per30, ...
             result.risetime, result.overshoot, result.settlingtime, ...
             result.divergence, result.SSE];

all_pass = true;
for i = 1:length(labels)
    e = exp_vals(i); c = comp_vals(i);
    if abs(e-c) < tol
        status = 'PASS';
    else
        status = sprintf('FAIL (%.2e)', e-c);
        all_pass = false;
    end
    fprintf('%-22s  %10.4f  %10.4f  %s\n', labels{i}, e, c, status);
end
fprintf('%s\n', repmat('-',58,1));
if all_pass
    fprintf('  All metrics matched within tolerance (%.0e).\n', tol);
else
    fprintf('  One or more metrics did NOT match.\n');
end
fprintf('=================================================================\n');

%% -----------------------------------------------------------------------
%  6. FOUR-GROUP SETTLING DIAGNOSTIC
%% -----------------------------------------------------------------------
fprintf('\n--- Settling group summary (threshold = 10.00 min) ---\n');
fprintf('  Group  Start (local)  End (local)  Duration   Result\n');
fprintf('  %-5d  %-13d  %-11d  %6.2f min  %s\n', ...
    1, 201,  334,  1.34, 'FAIL (trivially short)');
fprintf('  %-5d  %-13d  %-11d  %6.2f min  %s\n', ...
    2, 656,  1526, 8.70, 'FAIL (only 1.30 min short of threshold)');
fprintf('  %-5d  %-13d  %-11d  %6.2f min  %s\n', ...
    3, 1826, 3388, 15.62,'PASS → settle = 28.25 min');
fprintf('  %-5d  %-13d  %-11d  %6.2f min  %s\n', ...
    4, 3726, 8023, 42.97,'PASS (subsequent, does not change result)');
fprintf('\n');

%% -----------------------------------------------------------------------
%  7. HEMORRHAGE AND SSE CASE 2 DIAGNOSTIC
%% -----------------------------------------------------------------------
fprintf('--- Hemorrhage and SSE Case 2 summary ---\n');
fprintf('  Hem1: t = 24–26 min  (global idx 2401–2601)\n');
fprintf('  Hem2: t = 42–45 min  (global idx 4201–4501)\n');
fprintf('  Hem3: t = 88–92 min  (global idx 8801–9201)  ← last episode\n');
fprintf('  Post-Hem3 to CL end: ~8 min < 10 min → SSE Case 2 triggered\n');
fprintf('  While loop walks back through all 3 episodes (~6500 iterations)\n');
fprintf('  Converges to before Hem1 onset. Window: Phase 3, MAP=73.\n');
fprintf('  SSE = 73–70 = %.2f mmHg (correct regardless of convergence point\n', ...
        result.SSE);
fprintf('  since ALL pre-hemorrhage steady states are MAP=73 mmHg)\n\n');

%% -----------------------------------------------------------------------
%  8. PE INDEXING CHALLENGE DIAGNOSTIC
%% -----------------------------------------------------------------------
fprintf('--- PE indexing challenge ---\n');
PE_buggy   = (MAP(1 : n_eng) - target) / target * 100;
MDPE_buggy = median(PE_buggy);
fprintf('  Correct  MDPE (MAP(1001:10001)): %7.3f %%\n', result.MDPE);
fprintf('  Buggy    MDPE (MAP(1:9001)):      %7.3f %%\n', MDPE_buggy);
fprintf('  Pre-CL MAP=42 → PE=−40%% — clearly shifts buggy result\n\n');

%% -----------------------------------------------------------------------
%  9. WOBBLE DIAGNOSTIC
%% -----------------------------------------------------------------------
fprintf('--- Wobble diagnostic ---\n');
pct_at_73 = sum(MAP(ind_CL_start:ind_CL_end) == 73) / n_eng * 100;
fprintf('  Samples at MAP=73 exactly: %.1f%%  (< 50%%)\n', pct_at_73);
fprintf('  Phase 7 sinusoidal spreads |PE−MDPE| away from zero\n');
fprintf('  Wobble = %.4f %%  (non-zero confirmed)\n\n', result.Wobble);

%% -----------------------------------------------------------------------
%  10. SIGNAL PLOT
%% -----------------------------------------------------------------------
figure; hold on; box on;

% Shade pre-CL and CL regions
patch([0 10 10 0],       [38 38 92 92],[0.92 0.92 0.92],...
      'EdgeColor','none','DisplayName','Pre-CL (excluded)');
patch([10 100 100 10],   [38 38 92 92],[0.95 0.95 1.00],...
      'EdgeColor','none','DisplayName','CL engagement');

% Hemorrhage windows
patch([24 26 26 24],[38 38 92 92],[1.0 0.82 0.82],...
      'EdgeColor','none','DisplayName','Hem1 (t=24–26)');
patch([42 45 45 42],[38 38 92 92],[1.0 0.72 0.72],...
      'EdgeColor','none','DisplayName','Hem2 (t=42–45)');
patch([88 92 92 88],[38 38 92 92],[1.0 0.62 0.62],...
      'EdgeColor','none','DisplayName','Hem3 (t=88–92)');

% MAP and target
plot(Times, MAP, 'b-', 'LineWidth',1.5, 'DisplayName','MAP signal');
yline(target,     'k--','LineWidth',1.2,'DisplayName','Target (70 mmHg)');
yline(target*1.1, 'r:', 'LineWidth',1.0,'HandleVisibility','off');
yline(target*0.9, 'r:', 'LineWidth',1.0,'DisplayName','±10% band');

% Reference lines
xline(10, 'k-.','LineWidth',1.0,'DisplayName','CL start');
xline(Times(idx_T10_exp+ind_CL_start-1),'m--','LineWidth',1.0,'HandleVisibility','off');
xline(Times(idx_T90_exp+ind_CL_start-1),'m--','LineWidth',1.0,'HandleVisibility','off');
xline(result.settlingtime,'g--','LineWidth',1.2,...
      'DisplayName',sprintf('Settle=%.2f min (Group 3)',result.settlingtime));

% T10/T90 labels
text(Times(idx_T10_exp+ind_CL_start-1)+0.5,T10_lvl+4,'T_{10}','Color','m','FontSize',8);
text(Times(idx_T90_exp+ind_CL_start-1)+0.5,T90_lvl+4,'T_{90}','Color','m','FontSize',8);

% Group duration annotations
text(12, 88, 'G1:1.3m✗','Color',[0.5 0.5 0.5],'FontSize',7);
text(18, 88, 'G2:8.7m✗','Color',[0.8 0.4 0.0],'FontSize',7);
text(37, 88, 'G3:15.6m✓','Color',[0.0 0.6 0.0],'FontSize',7);
text(60, 88, 'G4:43.0m✓','Color',[0.0 0.4 0.8],'FontSize',7);

% Metric summary
text(11,83,sprintf('Overshoot = %.1f%%', result.overshoot),'Color','r','FontSize',8);
text(11,79,sprintf('MDPE = %.3f%%',      result.MDPE),      'Color','b','FontSize',8);
text(11,75,sprintf('Wobble = %.3f%%',    result.Wobble),     'Color','b','FontSize',8);
text(11,71,sprintf('SSE = %.2f mmHg (Case 2)',result.SSE),   'Color','b','FontSize',8);
text(3, 60, 'Pre-CL','Color',[0.5 0.5 0.5],'FontSize',8,'FontAngle','italic');

xlabel('Time (min)'); ylabel('MAP (mmHg)');
title('Validation Script 3 — Maximum Challenge — evaluate\_CL\_metrics');
legend('Location','southeast','FontSize',7,'NumColumns',2);
ylim([38 92]); xlim([0 100]); grid on;