% =========================================================
%  SCRIPT 1: Pipe Network Flow Analysis — Hardy-Cross Method
%  Subject:  Fluid Mechanics
%
%  Network Layout (7 pipes, 5 nodes, 2 loops):
%
%       [1]----P1----[2]----P2----[3]
%        |           |            |
%        P6          P4           P3
%        |           |            |
%       [5]----P5----[4]----------+
%
%  External Demands:
%    Node 1: Source (inflow Q_in)
%    Node 3: Demand 1
%    Node 5: Demand 2
%    Node 4: Demand 3
%
%  Method: Hardy-Cross iterative correction
%          Colebrook equation for friction factor (iterative)
%          Darcy-Weisbach for head loss
% =========================================================

clc; clear; close all;

fprintf('=======================================================\n');
fprintf('  SCRIPT 1: Hardy-Cross Pipe Network Solver\n');
fprintf('  Colebrook Friction Factors | Darcy-Weisbach\n');
fprintf('=======================================================\n\n');

% ---- Fluid Properties (Water at 20°C) ----
rho = 998;          % Density (kg/m³)
mu  = 1.002e-3;     % Dynamic viscosity (Pa·s)
nu  = mu / rho;     % Kinematic viscosity (m²/s)
g   = 9.81;         % Gravity (m/s²)

% ---- Pipe Properties ----
% Columns: [Length(m), Diameter(m), Roughness(mm)]
pipes = [
    500,  0.25,  0.15;   % Pipe 1: Node1-Node2
    400,  0.20,  0.15;   % Pipe 2: Node2-Node3
    350,  0.15,  0.26;   % Pipe 3: Node3-Node4
    300,  0.20,  0.15;   % Pipe 4: Node2-Node4
    450,  0.20,  0.15;   % Pipe 5: Node4-Node5
    380,  0.18,  0.26;   % Pipe 6: Node1-Node5
    320,  0.15,  0.26;   % Pipe 7: Node2-Node4 (parallel branch)
];

L  = pipes(:,1);
D  = pipes(:,2);
e  = pipes(:,3) * 1e-3;   % Convert mm to m
A  = pi/4 * D.^2;         % Cross-sectional areas

n_pipes = size(pipes,1);

% ---- Loop Structure ----
% Loop 1: Pipes 1, 4, 5, 6  (signs: +1 +4 -5 -6 based on clockwise positive)
% Loop 2: Pipes 2, 3, 4, 7
loop_pipes = {[1, 4, 5, 6], [2, 3, 4, 7]};
loop_signs = {[1, 1, -1, -1], [1, 1, -1, 1]};

% ---- Initial Flow Guesses (m³/s) ----
Q = [0.150; 0.080; 0.050; 0.070; 0.090; 0.060; 0.030];

% ---- Colebrook Solver (returns f given Re and e/D) ----
colebrook = @(Re, eD) fsolve(@(f) ...
    1/sqrt(f) + 2*log10(eD/3.7 + 2.51/(Re*sqrt(f))), ...
    0.02, optimset('Display','off'));

% ---- Hardy-Cross Iteration ----
max_iter  = 100;
tol       = 1e-6;   % m³/s convergence tolerance
iter      = 0;
dQ_hist   = zeros(max_iter, length(loop_pipes));

fprintf('Iterating Hardy-Cross...\n\n');
fprintf('%-8s %-15s %-15s\n', 'Iter', 'dQ_Loop1 (m³/s)', 'dQ_Loop2 (m³/s)');
fprintf('%s\n', repmat('-',1,42));

while iter < max_iter
    iter = iter + 1;
    converged = true;

    for loop = 1:length(loop_pipes)
        pids  = loop_pipes{loop};
        signs = loop_signs{loop};

        hf    = zeros(length(pids),1);
        dhfdQ = zeros(length(pids),1);

        for k = 1:length(pids)
            p  = pids(k);
            s  = signs(k);
            Qp = Q(p);
            V  = abs(Qp) / A(p);

            Re = V * D(p) / nu;
            if Re < 1
                Re = 1;  % avoid log(0)
            end

            eD = e(p) / D(p);

            % Friction factor via Colebrook (turbulent) or Hagen-Poiseuille
            if Re < 2300
                f = 64 / Re;
            elseif Re < 4000
                % Transition: interpolate
                f_lam  = 64 / 2300;
                f_turb = colebrook(4000, eD);
                f = f_lam + (Re-2300)/(4000-2300) * (f_turb - f_lam);
            else
                f = colebrook(Re, eD);
            end

            % Darcy-Weisbach head loss
            hf_p = f * L(p) / D(p) * V^2 / (2*g);
            hf(k)    = s * sign(Qp) * hf_p;
            dhfdQ(k) = 2 * f * L(p) / (D(p) * A(p)^2 * 2*g) * abs(Qp);
        end

        % Hardy-Cross correction
        dQ = -sum(hf) / sum(dhfdQ);
        dQ_hist(iter, loop) = dQ;

        % Apply correction
        for k = 1:length(pids)
            Q(pids(k)) = Q(pids(k)) + signs(k) * dQ;
        end

        if abs(dQ) > tol
            converged = false;
        end
    end

    if mod(iter,5)==0 || iter==1
        fprintf('%-8d %-15.6f %-15.6f\n', iter, ...
                dQ_hist(iter,1), dQ_hist(iter,2));
    end

    if converged
        break;
    end
end

fprintf('\nConverged in %d iterations.\n\n', iter);

% ---- Post-Processing ----
V_final  = abs(Q) ./ A;
Re_final = V_final .* D / nu;

% Friction factors at final solution
f_final = zeros(n_pipes,1);
for p = 1:n_pipes
    eD = e(p)/D(p);
    Re = Re_final(p);
    if Re < 2300
        f_final(p) = 64/Re;
    else
        f_final(p) = colebrook(max(Re,4000), eD);
    end
end

hf_final = f_final .* L ./ D .* V_final.^2 / (2*g);

% Flow regime classification
regime = cell(n_pipes,1);
for p = 1:n_pipes
    if Re_final(p) < 2300
        regime{p} = 'Laminar';
    elseif Re_final(p) < 4000
        regime{p} = 'Transitional';
    else
        regime{p} = 'Turbulent';
    end
end

% ---- Print Results Table ----
fprintf('%-6s %-10s %-10s %-10s %-10s %-10s %-14s\n', ...
        'Pipe', 'Q (L/s)', 'V (m/s)', 'Re', 'f', 'hf (m)', 'Regime');
fprintf('%s\n', repmat('-',1,74));
for p = 1:n_pipes
    fprintf('%-6d %-10.2f %-10.3f %-10.0f %-10.4f %-10.3f %-14s\n', ...
            p, Q(p)*1000, V_final(p), Re_final(p), f_final(p), hf_final(p), regime{p});
end

% ---- Save for other scripts ----
save('network_results.mat', 'Q','V_final','Re_final','f_final','hf_final', ...
     'regime','D','L','e','A','rho','mu','nu','g','pipes','n_pipes', ...
     'dQ_hist','iter','tol');

% ============================
%  PLOTS
% ============================

pipe_labels = {'P1','P2','P3','P4','P5','P6','P7'};
colors_regime = zeros(n_pipes,3);
for p=1:n_pipes
    if Re_final(p)<2300,      colors_regime(p,:)=[0.2 0.6 1.0];   % blue=laminar
    elseif Re_final(p)<4000,  colors_regime(p,:)=[1.0 0.7 0.0];   % orange=trans
    else,                     colors_regime(p,:)=[0.9 0.2 0.2];   % red=turbulent
    end
end

% --- Plot 1: Flow Rates ---
figure('Name','Network Flow Rates','Color','white','Position',[50 400 700 420]);
b = bar(Q*1000, 'FaceColor','flat');
b.CData = colors_regime;
hold on;
yline(0,'k--','LineWidth',1);
set(gca,'XTickLabel',pipe_labels,'FontSize',12);
xlabel('Pipe','FontSize',13,'FontWeight','bold');
ylabel('Flow Rate  Q  (L/s)','FontSize',13,'FontWeight','bold');
title({'Hardy-Cross Converged Flow Rates'; ...
       'Blue=Laminar | Orange=Transitional | Red=Turbulent'}, ...
      'FontSize',13,'FontWeight','bold');
grid on; grid minor;
for p=1:n_pipes
    text(p, Q(p)*1000 + 0.5*sign(Q(p)), sprintf('%.1f',Q(p)*1000), ...
         'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
end
box on;

% --- Plot 2: Reynolds Number with regime bands ---
figure('Name','Reynolds Numbers','Color','white','Position',[100 350 750 450]);
bar(Re_final, 'FaceColor','flat','CData',colors_regime);
hold on;
yline(2300, 'b--', 'LineWidth',2);
yline(4000, 'r--', 'LineWidth',2);
text(0.6, 2500, 'Re = 2300  (Laminar limit)', 'FontSize',10,'Color','b');
text(0.6, 4300, 'Re = 4000  (Turbulent onset)', 'FontSize',10,'Color','r');
set(gca,'XTickLabel',pipe_labels,'FontSize',12);
xlabel('Pipe','FontSize',13,'FontWeight','bold');
ylabel('Reynolds Number  Re','FontSize',13,'FontWeight','bold');
title('Reynolds Number per Pipe — Flow Regime Classification','FontSize',13,'FontWeight','bold');
grid on; grid minor; box on;

% --- Plot 3: Head Loss ---
figure('Name','Head Losses','Color','white','Position',[150 300 700 420]);
bar(hf_final,'FaceColor',[0.85 0.33 0.10]);
set(gca,'XTickLabel',pipe_labels,'FontSize',12);
xlabel('Pipe','FontSize',13,'FontWeight','bold');
ylabel('Head Loss  h_f  (m)','FontSize',13,'FontWeight','bold');
title('Darcy-Weisbach Head Losses Across Network','FontSize',13,'FontWeight','bold');
grid on; grid minor; box on;
for p=1:n_pipes
    text(p, hf_final(p)+0.05, sprintf('%.2f m',hf_final(p)), ...
         'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color','k');
end

% --- Plot 4: Hardy-Cross Convergence ---
figure('Name','Convergence','Color','white','Position',[200 250 700 420]);
semilogy(1:iter, abs(dQ_hist(1:iter,1))*1000, 'b-o', 'LineWidth',2,'MarkerSize',5); hold on;
semilogy(1:iter, abs(dQ_hist(1:iter,2))*1000, 'r-s', 'LineWidth',2,'MarkerSize',5);
yline(tol*1000,'k--','LineWidth',1.5);
text(iter*0.6, tol*1000*3, sprintf('Tolerance = %.2e L/s', tol*1000), ...
     'FontSize',10,'Color','k');
xlabel('Iteration','FontSize',13,'FontWeight','bold');
ylabel('|\DeltaQ|  (L/s)  — log scale','FontSize',13,'FontWeight','bold');
title('Hardy-Cross Convergence History','FontSize',13,'FontWeight','bold');
legend('Loop 1','Loop 2','Tolerance','Location','northeast','FontSize',11);
grid on; box on;

fprintf('\n[Script 1 complete] Network solved. Results saved to network_results.mat\n');
