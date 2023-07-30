%% LOAD DATA

data.MultiRTAC_ms = dlmread("stats_MultiRTAC__ms.csv", ',', 1, 0);
data.MultiRTAC_ss = dlmread("stats_MultiRTAC__ss.csv", ',', 1, 0);

data.QuadcopterFull_ms = dlmread("stats_QuadcopterFull__ms.csv", ',', 1, 0);
data.QuadcopterFull_ss = dlmread("stats_QuadcopterFull__ss.csv", ',', 1, 0);

data.NonlinearOCP1b_ms = dlmread("stats_NonlinearOCP1b__ms.csv", ',', 1, 0);
data.NonlinearOCP1b_ss = dlmread("stats_NonlinearOCP1b__ss.csv", ',', 1, 0);

%% CALCULATE SPEED-UP

% MultiRTAC:
multi_rtac_time_iter_ms = data.MultiRTAC_ms(2,9)/data.MultiRTAC_ms(2,5);
multi_rtac_time_iter_ss = data.MultiRTAC_ss(1,8)/data.MultiRTAC_ss(1,4);
multi_rtac_speedup = multi_rtac_time_iter_ss/multi_rtac_time_iter_ms;

% Quadcopter Full:
quadcopter_full_time_iter_ms = data.QuadcopterFull_ms(2,9)/data.QuadcopterFull_ms(2,5);
quadcopter_full_time_iter_ss = data.QuadcopterFull_ss(2,8);
quadcopter_full_speedup = quadcopter_full_time_iter_ss/quadcopter_full_time_iter_ms;

% NonlinearOCP1b:
nonlinear_time_iter_ms = data.NonlinearOCP1b_ms(3,9)/data.NonlinearOCP1b_ms(3,5);
nonlinear_time_iter_ss = data.NonlinearOCP1b_ss(2,8)/data.MultiRTAC_ss(1,4);
nonlinear_speedup = nonlinear_time_iter_ss/nonlinear_time_iter_ms;


