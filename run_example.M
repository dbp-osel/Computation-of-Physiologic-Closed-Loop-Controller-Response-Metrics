clear all
clc
close all
load('test_data.mat')
PCLC_char = evaluate_CL_metrics(ind_CL_start,ind_CL_end,MAP_Values,MAP_Times,target,Hemorrhage_Values)


figure;
% --- Top plot: MAP vs Time ---
subplot(2,1,1)
plot(MAP_Times, MAP_Values, 'b', 'LineWidth', 1.5);
hold on;
% Target line
yline(target, 'r--', 'LineWidth', 2);
% CL interval markers
xline(MAP_Times(ind_CL_start), 'g--', 'PCLC Start', 'LineWidth', 1.5);
xline(MAP_Times(ind_CL_end),   'm--', 'PCLC End',   'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('MAP (mmHg)');
%title('MAP vs Time');
legend('MAP', 'Target');
grid on;
% --- Bottom plot: Hemorrhage presence ---
subplot(2,1,2)
stairs(MAP_Times, Hemorrhage_Values, 'k', 'LineWidth', 1.5);
hold on;
% CL interval markers
xline(MAP_Times(ind_CL_start), 'g--', 'PCLC Start', 'LineWidth', 1.5);
xline(MAP_Times(ind_CL_end),   'm--', 'PCLC End',   'LineWidth', 1.5);
xlabel('Time');
ylabel('Hemorrhage (l/min)');
%title('Hemorrhage Presence');
grid on;