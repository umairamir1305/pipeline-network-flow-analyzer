% =========================================================
%  MASTER RUNNER Pipeline Network Flow Analysis Suite
%
%  Author : Umair
%  Project: Pipeline Network Flow Analyzer with
%           Pump Selection & Cavitation Risk Assessment
%
%  Run this file. All scripts execute in order.
%  Requires: Optimization Toolbox (for fsolve/fzero)
% =========================================================

clc; clear; close all;

fprintf('\n');
fprintf('#######################################################\n');
fprintf('#                                                     #\n');
fprintf('#   PIPELINE NETWORK FLOW ANALYSIS SUITE             #\n');
fprintf('#   Hardy-Cross | Moody | Pump Curves | Cavitation   #\n');
fprintf('#   Boundary Layer | Parametric Study                #\n');
fprintf('#                                                     #\n');
fprintf('#######################################################\n\n');

% Check toolbox availability
if ~license('test','optimization_toolbox')
    error(['Optimization Toolbox required for fsolve/fzero.\n' ...
           'Available via NUST MATLAB campus license.']);
end

fprintf('All dependencies OK. Starting analysis...\n\n');
fprintf('-------------------------------------------------------\n');

%% Script 1
fprintf('\n>>> SCRIPT 1: Hardy-Cross Network Solver\n');
fprintf('    Colebrook friction factors | 7 pipes | 2 loops\n\n');
run('script1_hardy_cross_network.m');
fprintf('-------------------------------------------------------\n');

%% Script 2
fprintf('\n>>> SCRIPT 2: Moody Diagram\n');
fprintf('    Full diagram + network operating points overlaid\n\n');
run('script2_moody_diagram.m');
fprintf('-------------------------------------------------------\n');

%% Script 3
fprintf('\n>>> SCRIPT 3: Pump Selection & Cavitation\n');
fprintf('    3 pump curves | NPSH analysis | cavitation flags\n\n');
run('script3_pump_cavitation.m');
fprintf('-------------------------------------------------------\n');

%% Script 4
fprintf('\n>>> SCRIPT 4: Boundary Layer + Parametric Study\n');
fprintf('    Entry length | velocity profiles | diameter sweep\n\n');
run('script4_boundary_layer_parametric.m');
fprintf('-------------------------------------------------------\n');

fprintf('\n#######################################################\n');
fprintf('#   ALL SCRIPTS COMPLETE                              #\n');
fprintf('#   Total figures generated: ~12                     #\n');
fprintf('#   Intermediate data: network_results.mat           #\n');
fprintf('#######################################################\n\n');
