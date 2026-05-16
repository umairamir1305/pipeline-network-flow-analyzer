% =========================================================
%  SCRIPT 3: Pump Selection, Operating Point & Cavitation
%  Subject:  Turbomachinery + Bernoulli + Fluid Mechanics
%
%  Theory:
%    System curve:   H_sys = H_static + R * Q^2
%      where R = sum of (f*L/D * 1/(2g*A^2)) for series pipes
%
%    Pump curve:     H_pump = H_shutoff - a*Q - b*Q^2
%      (manufacturer polynomial fit)
%
%    Operating point: intersection of pump and system curves
%      solved via fzero
%
%    Cavitation (NPSH):
%      NPSH_available = (p_atm - p_vapour)/(rho*g) + V_s^2/(2g) - z_s
%      NPSH_required  = from pump curve (manufacturer data)
%      Cavitation occurs if NPSH_a < NPSH_r
% =========================================================

clc; close all;

fprintf('=======================================================\n');
fprintf('  SCRIPT 3: Pump Operating Point & Cavitation\n');
fprintf('=======================================================\n\n');

load('network_results.mat');

% ---- System Parameters ----
H_static   = 15;       % Static head (m) — elevation difference source to delivery
p_atm      = 101325;   % Atmospheric pressure (Pa)
p_vapour   = 2337;     % Vapour pressure of water at 20°C (Pa)
z_s        = -2.0;     % Pump inlet elevation relative to water surface (m)
                       % negative = pump is ABOVE the water surface (suction lift)

% ---- System Resistance (main supply branch: pipes 1+2) ----
% R = sum of f*L/(D * 2g * A^2) for pipes in series
% We compute it empirically from Script 1 results
% Use pipes 1 and 2 as the main branch
main_pipes = [1, 2];
Q_main     = mean(abs(Q(main_pipes)));     % approximate main flow

R_total = 0;
for p = main_pipes
    Vp  = abs(Q(p)) / A(p);
    Re  = Vp * D(p) / nu;
    eD  = e(p) / D(p);
    if Re < 2300
        fp = 64/Re;
    else
        fp = f_final(p);
    end
    R_total = R_total + fp * L(p) / D(p) / (2*g*A(p)^2);
end

% System curve: H_sys(Q) = H_static + R_total * Q^2
Q_sys  = linspace(0, 0.35, 300);     % m³/s
H_sys  = H_static + R_total * Q_sys.^2;

% ---- Pump Curve (Centrifugal Pump) ----
% H = H0 - a*Q - b*Q^2
% Three pumps to compare
pumps(1).name  = 'Pump A (Standard)';
pumps(1).H0    = 45;   pumps(1).a = 20;   pumps(1).b = 350;
pumps(1).eta0  = 0.72; pumps(1).NPSH_r = 3.5;   % required NPSH (m)
pumps(1).color = [0.2 0.5 0.9];

pumps(2).name  = 'Pump B (High-Head)';
pumps(2).H0    = 60;   pumps(2).a = 30;   pumps(2).b = 500;
pumps(2).eta0  = 0.68; pumps(2).NPSH_r = 5.0;
pumps(2).color = [0.9 0.3 0.2];

pumps(3).name  = 'Pump C (High-Flow)';
pumps(3).H0    = 35;   pumps(3).a = 10;   pumps(3).b = 200;
pumps(3).eta0  = 0.75; pumps(3).NPSH_r = 2.8;
pumps(3).color = [0.1 0.75 0.4];

fprintf('%-25s %-14s %-12s %-12s %-12s %-10s\n', ...
        'Pump','Q_op (L/s)','H_op (m)','Power (kW)','Efficiency','Cavitation');
fprintf('%s\n',repmat('-',1,90));

% ---- Figure: Pump & System Curves ----
figure('Name','Pump Curves & System Curve','Color','white','Position',[50 50 850 550]);
hold on;

% System curve
plot(Q_sys*1000, H_sys, 'k-', 'LineWidth', 3);
text(Q_sys(end)*1000*0.92, H_sys(end)+1.5, 'System Curve', ...
     'FontSize',11,'FontWeight','bold','Color','k');

op_points = zeros(length(pumps), 2);

for i = 1:length(pumps)
    pm   = pumps(i);
    H_pump_fn = @(Q) pm.H0 - pm.a*Q - pm.b*Q.^2;
    H_pump    = H_pump_fn(Q_sys);
    H_pump(H_pump < 0) = NaN;

    % Plot pump curve
    plot(Q_sys*1000, H_pump, '-', 'Color', pm.color, 'LineWidth', 2.5);

    % Shut-off head label
    text(2, pm.H0 + 0.5, sprintf('H_0=%.0fm', pm.H0), ...
         'FontSize', 9, 'Color', pm.color);

    % Find operating point: H_pump(Q) = H_sys(Q)
    f_op = @(Q) (pm.H0 - pm.a*Q - pm.b*Q^2) - (H_static + R_total*Q^2);
    try
        Q_op = fzero(f_op, [0.01, 0.30]);
        H_op = H_static + R_total * Q_op^2;
    catch
        Q_op = NaN; H_op = NaN;
    end

    op_points(i,:) = [Q_op, H_op];

    % Mark operating point
    if ~isnan(Q_op)
        plot(Q_op*1000, H_op, 'o', 'MarkerSize', 14, ...
             'MarkerFaceColor', pm.color, 'MarkerEdgeColor','k','LineWidth',1.5);
        text(Q_op*1000 + 2, H_op + 1.5, ...
             sprintf('OP: %.1f L/s\n    %.1f m', Q_op*1000, H_op), ...
             'FontSize', 9.5, 'Color', pm.color, 'FontWeight','bold');
    end

    % Power and efficiency
    V_s      = Q_op / (pi/4 * 0.10^2);     % suction pipe velocity (D=100mm assumed)
    NPSH_a   = (p_atm - p_vapour)/(rho*g) + V_s^2/(2*g) - abs(z_s);
    Power_kW = rho * g * Q_op * H_op / pm.eta0 / 1000;

    if NPSH_a >= pm.NPSH_r
        cav_str = 'SAFE';
    else
        cav_str = '*** CAVITATION ***';
    end

    fprintf('%-25s %-14.2f %-12.2f %-12.2f %-12.2f %-10s\n', ...
            pm.name, Q_op*1000, H_op, Power_kW, pm.eta0*100, cav_str);
end

xlabel('Flow Rate,  Q  (L/s)',  'FontSize',13,'FontWeight','bold');
ylabel('Head,  H  (m)',         'FontSize',13,'FontWeight','bold');
title({'Pump Curves vs System Curve — Operating Point Determination'; ...
       sprintf('H_{static} = %.0f m  |  System resistance R = %.1f m/(m^3/s)^2', ...
               H_static, R_total)}, 'FontSize',13,'FontWeight','bold');
legend('System curve', pumps(1).name, pumps(2).name, pumps(3).name, ...
       'Location','northeast','FontSize',10);
xlim([0, 320]); ylim([0, 75]);
grid on; grid minor; box on;
set(gca,'FontSize',11,'LineWidth',1.2);

% ---- Figure: NPSH & Cavitation Analysis ----
figure('Name','Cavitation Analysis','Color','white','Position',[100 100 800 500]);

Q_cav   = linspace(0.01, 0.30, 300);
V_s_arr = Q_cav / (pi/4 * 0.10^2);
NPSH_a  = (p_atm - p_vapour)/(rho*g) + V_s_arr.^2/(2*g) - abs(z_s);

hold on;
plot(Q_cav*1000, NPSH_a, 'b-', 'LineWidth', 3);

pump_colors = {pumps(1).color, pumps(2).color, pumps(3).color};
for i = 1:length(pumps)
    NPSH_r_line = pumps(i).NPSH_r * ones(size(Q_cav));
    plot(Q_cav*1000, NPSH_r_line, '--', 'Color', pump_colors{i}, 'LineWidth', 2);

    % Cavitation onset point
    f_cav = @(Q) (p_atm-p_vapour)/(rho*g) + (Q/(pi/4*0.10^2))^2/(2*g) ...
                 - abs(z_s) - pumps(i).NPSH_r;
    try
        Q_cav_onset = fzero(f_cav, [0.01, 0.30]);
        NPSH_onset  = pumps(i).NPSH_r;
        plot(Q_cav_onset*1000, NPSH_onset, 'x', 'MarkerSize', 14, ...
             'LineWidth', 3, 'Color', pump_colors{i});
        text(Q_cav_onset*1000+2, NPSH_onset+0.15, ...
             sprintf('Cavitation onset\n@ %.1f L/s', Q_cav_onset*1000), ...
             'FontSize',9,'Color',pump_colors{i},'FontWeight','bold');
    catch
    end

    % Mark operating point
    Q_op = op_points(i,1);
    if ~isnan(Q_op)
        V_s_op  = Q_op/(pi/4*0.10^2);
        NPSH_op = (p_atm-p_vapour)/(rho*g) + V_s_op^2/(2*g) - abs(z_s);
        plot(Q_op*1000, NPSH_op, 'o', 'MarkerSize',10, ...
             'MarkerFaceColor',pump_colors{i},'MarkerEdgeColor','k');
    end
end

% Shade cavitation danger zone
ylims = ylim;
patch([0 300 300 0],[ylims(1) ylims(1) 2.0 2.0], ...
      [1 0.8 0.8],'FaceAlpha',0.3,'EdgeColor','none');
text(10, 1.2,'Cavitation Risk Zone','FontSize',10,'Color',[0.8 0 0],'FontWeight','bold');

xlabel('Flow Rate,  Q  (L/s)',        'FontSize',13,'FontWeight','bold');
ylabel('NPSH  (m)',                    'FontSize',13,'FontWeight','bold');
title({'Net Positive Suction Head — Cavitation Assessment'; ...
       sprintf('Suction lift z_s = %.1f m  |  T = 20°C', abs(z_s))}, ...
      'FontSize',13,'FontWeight','bold');
legend('NPSH_{available}', ...
       [pumps(1).name ' NPSH_r'], [pumps(2).name ' NPSH_r'], [pumps(3).name ' NPSH_r'], ...
       'Location','southeast','FontSize',10);
grid on; grid minor; box on;
set(gca,'FontSize',11,'LineWidth',1.2);

fprintf('\n[Script 3 complete] Pump selection and cavitation analysis done.\n');
