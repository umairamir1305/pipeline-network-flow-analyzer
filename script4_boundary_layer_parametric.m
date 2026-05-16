% =========================================================
%  SCRIPT 4: Boundary Layer Development & Parametric Study
%  Subject:  Boundary Layer Theory + Full System
%
%  Part A: Boundary Layer development in pipe entry region
%    - Hydrodynamic entry length for laminar & turbulent flow
%    - Velocity profile evolution from flat to parabolic (laminar)
%      or turbulent power-law profile
%    - Wall shear stress distribution along entry length
%
%  Part B: Parametric Study
%    - Vary pipe diameter of main supply pipe
%    - Show effect on: Re, f, hf, flow distribution, pump demand
% =========================================================

clc; close all;

fprintf('=======================================================\n');
fprintf('  SCRIPT 4: Boundary Layer + Parametric Study\n');
fprintf('=======================================================\n\n');

load('network_results.mat');

% =====================================================
%  PART A: BOUNDARY LAYER / ENTRY LENGTH ANALYSIS
% =====================================================

fprintf('--- PART A: Entry Length & Velocity Profile ---\n\n');

% Analyse Pipe 1 (largest diameter, highest Re)
p_sel  = 1;
D_p    = D(p_sel);
R_p    = D_p / 2;
Q_p    = abs(Q(p_sel));
V_avg  = Q_p / A(p_sel);
Re_p   = Re_final(p_sel);
f_p    = f_final(p_sel);

fprintf('Analysing Pipe 1:\n');
fprintf('  D = %.0f mm | V_avg = %.3f m/s | Re = %.0f\n\n', D_p*1000, V_avg, Re_p);

% Entry lengths
if Re_p < 2300
    L_entry = 0.06 * Re_p * D_p;    % Laminar: L_h = 0.06*Re*D
    profile_type = 'Laminar (Parabolic)';
else
    L_entry = 4.4 * Re_p^(1/6) * D_p;  % Turbulent: L_h = 4.4*Re^(1/6)*D
    profile_type = 'Turbulent (Power-law 1/7)';
end

fprintf('  Entry length L_h = %.2f m  (%s)\n', L_entry, profile_type);
fprintf('  Pipe length  L   = %.0f m\n\n', L(p_sel));

if L_entry < L(p_sel)
    fprintf('  Flow is FULLY DEVELOPED for most of the pipe.\n\n');
else
    fprintf('  WARNING: Flow may NOT be fully developed!\n\n');
end

% ---- Velocity profiles at different x/L_entry stations ----
r_norm = linspace(0, 1, 200);   % r/R from 0 (centre) to 1 (wall)
x_stations = [0.01, 0.1, 0.3, 0.6, 1.0];  % x/L_entry

figure('Name','Velocity Profile Development','Color','white','Position',[50 400 750 550]);
hold on;

cmap2 = cool(length(x_stations));
legend_str = cell(length(x_stations),1);

for k = 1:length(x_stations)
    xL = x_stations(k);

    if Re_p < 2300
        % Laminar: profile develops from flat toward parabolic
        % Approximate: blend flat and parabolic by development fraction
        dev = min(xL, 1.0);
        V_flat     = V_avg * ones(size(r_norm));
        V_parabola = 2 * V_avg * (1 - r_norm.^2);
        V_profile  = (1-dev)*V_flat + dev*V_parabola;
    else
        % Turbulent: 1/7 power law, sharpens slightly
        % Near entry: flatter; developed: 1/7 power
        n  = 7 + 2.5*min(xL,1.0);     % n increases toward developed value
        V_profile = V_avg * (n+1)*(2*n+1)/(2*n^2) * (1 - r_norm).^(1/n);
        V_profile(end) = 0;           % no-slip at wall
    end

    plot(V_profile, r_norm, '-', 'Color', cmap2(k,:), 'LineWidth', 2.2);
    legend_str{k} = sprintf('x/L_h = %.2f', xL);
end

% Mark centreline and wall
xline_val = V_avg;
plot([xline_val xline_val],[0 1],'k--','LineWidth',1.5);
text(xline_val+0.02, 0.95, 'V_{avg}','FontSize',10,'Color','k');

xlabel('Velocity  u  (m/s)',     'FontSize',13,'FontWeight','bold');
ylabel('r/R  (0=centre, 1=wall)','FontSize',13,'FontWeight','bold');
title({sprintf('Velocity Profile Development — Pipe 1 (Re = %.0f)', Re_p); ...
       sprintf('%s | L_h = %.1f m', profile_type, L_entry)}, ...
      'FontSize',13,'FontWeight','bold');
legend(legend_str,'Location','southwest','FontSize',10);
grid on; grid minor; box on;
set(gca,'FontSize',11,'LineWidth',1.2);
set(gca,'YDir','normal');

% ---- Wall Shear Stress along entry length ----
x_arr = linspace(0.01, L(p_sel), 400);

if Re_p < 2300
    % Laminar developing: tau_w = 3.44/sqrt(x/D/Re) * (rho*V^2/2) + fully dev value
    x_star    = x_arr / (D_p * Re_p);    % dimensionless x
    Cf_dev    = 16 / Re_p;               % fully developed Cf = f/4
    Cf_entry  = 3.44 ./ sqrt(pi * x_star) + Cf_dev;
    Cf_entry(Cf_entry > 10*Cf_dev) = 10*Cf_dev;
    tau_w     = Cf_entry * (0.5 * rho * V_avg^2);
else
    % Turbulent: tau_w = f/8 * rho * V^2 (constant in fully developed region)
    % Entry region slightly higher
    x_star    = x_arr / D_p;
    Cf_fd     = f_p / 8;
    decay     = 1 + 3 * exp(-x_arr / L_entry);
    tau_w     = Cf_fd * decay * rho * V_avg^2;
end

figure('Name','Wall Shear Stress','Color','white','Position',[100 350 750 450]);
plot(x_arr, tau_w, 'r-','LineWidth',2.5); hold on;
xline(L_entry,'b--','LineWidth',2);
text(L_entry*1.02, max(tau_w)*0.85, sprintf('L_h = %.1f m', L_entry), ...
     'FontSize',10,'Color','b','FontWeight','bold');

xlabel('Distance from entry,  x  (m)',  'FontSize',13,'FontWeight','bold');
ylabel('Wall Shear Stress  \tau_w  (Pa)','FontSize',13,'FontWeight','bold');
title({'Wall Shear Stress Distribution Along Pipe Entry Region'; ...
       sprintf('Pipe 1 | D = %.0f mm | Re = %.0f', D_p*1000, Re_p)}, ...
      'FontSize',13,'FontWeight','bold');
grid on; grid minor; box on;
set(gca,'FontSize',11,'LineWidth',1.2);

% =====================================================
%  PART B: PARAMETRIC STUDY — Diameter vs System Response
% =====================================================

fprintf('\n--- PART B: Parametric Study — Pipe Diameter Effect ---\n\n');

D_range   = linspace(0.10, 0.40, 60);   % Vary main pipe diameter 100-400mm
hf_arr    = zeros(size(D_range));
Re_arr    = zeros(size(D_range));
f_arr     = zeros(size(D_range));
V_arr     = zeros(size(D_range));
Q_fixed   = abs(Q(1));   % Hold flow fixed at Script 1 solution

colebrook = @(Re, eD) fsolve(@(f) ...
    1/sqrt(f) + 2*log10(eD/3.7 + 2.51/(Re*sqrt(f))), ...
    0.02, optimset('Display','off'));

for i = 1:length(D_range)
    Di  = D_range(i);
    Ai  = pi/4 * Di^2;
    Vi  = Q_fixed / Ai;
    Rei = Vi * Di / nu;
    eDi = e(1) / Di;

    if Rei < 2300
        fi = 64/Rei;
    elseif Rei < 4000
        fi = 64/2300 + (Rei-2300)/1700 * (colebrook(4000,eDi) - 64/2300);
    else
        fi = colebrook(Rei, eDi);
    end

    hf_arr(i) = fi * L(1) / Di * Vi^2 / (2*g);
    Re_arr(i) = Rei;
    f_arr(i)  = fi;
    V_arr(i)  = Vi;
end

% Find diameter for Re = 2300 (laminar-turbulent transition)
D_trans = interp1(Re_arr, D_range, 2300, 'linear','extrap');
fprintf('Transition diameter (Re=2300) for Q=%.1f L/s: D = %.1f mm\n', ...
        Q_fixed*1000, D_trans*1000);

figure('Name','Parametric Study — Diameter','Color','white','Position',[200 100 900 600]);

subplot(2,2,1);
plot(D_range*1000, Re_arr,'b-','LineWidth',2.2); hold on;
yline(2300,'b--','LineWidth',1.5); yline(4000,'r--','LineWidth',1.5);
text(D_range(end)*1000*0.7, 2600,'Re=2300','FontSize',9,'Color','b');
text(D_range(end)*1000*0.7, 4300,'Re=4000','FontSize',9,'Color','r');
xlabel('Diameter (mm)','FontSize',11,'FontWeight','bold');
ylabel('Reynolds Number','FontSize',11,'FontWeight','bold');
title('Re vs Diameter','FontSize',12,'FontWeight','bold');
grid on; grid minor; box on;

subplot(2,2,2);
plot(D_range*1000, f_arr,'r-','LineWidth',2.2);
xlabel('Diameter (mm)','FontSize',11,'FontWeight','bold');
ylabel('Friction Factor f','FontSize',11,'FontWeight','bold');
title('Friction Factor vs Diameter','FontSize',12,'FontWeight','bold');
grid on; grid minor; box on;

subplot(2,2,3);
semilogy(D_range*1000, hf_arr,'m-','LineWidth',2.2); hold on;
plot(D(1)*1000, hf_final(1),'ko','MarkerFaceColor','m','MarkerSize',10);
text(D(1)*1000+3, hf_final(1)*1.3, sprintf('Current\nD=%.0fmm\nh_f=%.1fm', ...
     D(1)*1000, hf_final(1)),'FontSize',9,'Color','m','FontWeight','bold');
xlabel('Diameter (mm)','FontSize',11,'FontWeight','bold');
ylabel('Head Loss h_f (m) — log','FontSize',11,'FontWeight','bold');
title('Head Loss vs Diameter','FontSize',12,'FontWeight','bold');
grid on; grid minor; box on;

subplot(2,2,4);
plot(D_range*1000, V_arr,'g-','LineWidth',2.2); hold on;
yline(3.0,'r--','LineWidth',1.5);  % typical max velocity for water mains
text(D_range(3)*1000, 3.15,'V_{max,recommended}=3 m/s','FontSize',9,'Color','r');
xlabel('Diameter (mm)','FontSize',11,'FontWeight','bold');
ylabel('Velocity (m/s)','FontSize',11,'FontWeight','bold');
title('Velocity vs Diameter','FontSize',12,'FontWeight','bold');
grid on; grid minor; box on;

sgtitle({'Parametric Study: Effect of Main Pipe Diameter on Flow Characteristics'; ...
         sprintf('Q = %.1f L/s fixed | Pipe 1 | L = %.0f m | e = %.2f mm', ...
                 Q_fixed*1000, L(1), e(1)*1000)}, ...
        'FontSize',13,'FontWeight','bold');

fprintf('\n[Script 4 complete] Boundary layer and parametric study done.\n');
