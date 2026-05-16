% =========================================================
%  SCRIPT 2: Moody Diagram Generator + Pipe Operating Points
%  Subject:  Viscous Flow / Reynolds Number
%
%  Generates a full Moody diagram and overlays the actual
%  operating points of every pipe in the network.
%
%  Friction factor correlations used:
%    Laminar:     f = 64/Re
%    Turbulent:   Colebrook-White equation (implicit, solved iteratively)
%    Smooth pipe: Prandtl smooth-pipe law
% =========================================================

clc; close all;

fprintf('=======================================================\n');
fprintf('  SCRIPT 2: Moody Diagram + Network Operating Points\n');
fprintf('=======================================================\n\n');

% ---- Load network results ----
if ~exist('network_results.mat','file')
    error('Run Script 1 first to generate network_results.mat');
end
load('network_results.mat');

% ---- Colebrook solver ----
colebrook = @(Re, eD) fsolve(@(f) ...
    1/sqrt(f) + 2*log10(eD/3.7 + 2.51/(Re*sqrt(f))), ...
    0.02, optimset('Display','off'));

% ---- Build Moody Diagram ----
Re_range   = logspace(2.7, 8, 600);
eD_values  = [0, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 5e-3, 1e-2];
eD_labels  = {'Smooth','1e-6','5e-6','1e-5','5e-5','1e-4', ...
               '2e-4','5e-4','1e-3','5e-3','1e-2'};

% Color map for relative roughness lines
cmap = parula(length(eD_values));

figure('Name','Moody Diagram','Color','white','Position',[50 50 950 650]);
hold on;

% Laminar line
Re_lam = logspace(2.7, log10(2300), 100);
f_lam  = 64 ./ Re_lam;
loglog(Re_lam, f_lam, 'k-', 'LineWidth', 3);

% Critical zone shading
patch([2300 4000 4000 2300],[0.005 0.005 0.1 0.1], ...
      [0.9 0.9 0.9],'FaceAlpha',0.5,'EdgeColor','none');
text(3000, 0.006, {'Transition','Zone'}, 'FontSize',8, ...
     'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

% Turbulent lines for each relative roughness
for i = 1:length(eD_values)
    eD  = eD_values(i);
    f_t = zeros(size(Re_range));
    for j = 1:length(Re_range)
        Re = Re_range(j);
        if Re < 2300
            continue;
        elseif eD == 0
            % Smooth pipe — Prandtl
            f_t(j) = fsolve(@(f) 1/sqrt(f) - 2*log10(Re*sqrt(f)) + 0.8, ...
                             0.02, optimset('Display','off'));
        else
            f_t(j) = colebrook(Re, eD);
        end
    end
    idx = Re_range >= 4000;
    loglog(Re_range(idx), f_t(idx), '-', 'Color', cmap(i,:), 'LineWidth', 1.4);

    % Label at right edge
    f_end = f_t(end);
    if f_end > 0
        text(Re_range(end)*1.05, f_end, eD_labels{i}, ...
             'FontSize', 7.5, 'Color', cmap(i,:), 'FontWeight','bold');
    end
end

% ---- Overlay Network Pipe Operating Points ----
pipe_labels = {'P1','P2','P3','P4','P5','P6','P7'};
markers = {'o','s','^','d','v','p','h'};

for p = 1:n_pipes
    Re_p = Re_final(p);
    f_p  = f_final(p);
    loglog(Re_p, f_p, markers{mod(p-1,7)+1}, ...
           'MarkerSize', 12, 'LineWidth', 2, ...
           'MarkerFaceColor', [0.95 0.3 0.1], ...
           'MarkerEdgeColor', 'k');
    text(Re_p * 1.08, f_p, pipe_labels{p}, ...
         'FontSize', 10, 'FontWeight','bold', 'Color',[0.8 0.1 0.1]);
end

% ---- Annotations ----
text(500,   0.055, 'Laminar', 'FontSize',11,'FontWeight','bold','Color','k');
text(1.5e7, 0.055, 'e/D values \rightarrow', 'FontSize',9,'Color',[0.3 0.3 0.3]);

xlabel('Reynolds Number,  Re', 'FontSize',13,'FontWeight','bold');
ylabel('Darcy Friction Factor,  f', 'FontSize',13,'FontWeight','bold');
title({'Moody Diagram — Friction Factor vs Reynolds Number'; ...
       'Network pipe operating points overlaid (red markers)'}, ...
      'FontSize',13,'FontWeight','bold');

xlim([5e2, 2e8]);
ylim([0.005, 0.1]);
set(gca,'FontSize',11,'LineWidth',1.2,'XGrid','on','YGrid','on', ...
        'GridAlpha',0.25,'MinorGridAlpha',0.1);
box on;

% ---- Print operating summary ----
fprintf('%-6s %-12s %-10s %-10s %-12s\n','Pipe','Re','f','e/D','Regime');
fprintf('%s\n',repmat('-',1,54));
eD_net = pipes(:,3)*1e-3 ./ pipes(:,2);
for p=1:n_pipes
    if Re_final(p)<2300,      reg='Laminar';
    elseif Re_final(p)<4000,  reg='Transitional';
    else,                     reg='Turbulent';
    end
    fprintf('%-6s %-12.0f %-10.4f %-10.5f %-12s\n', ...
            pipe_labels{p}, Re_final(p), f_final(p), eD_net(p), reg);
end

fprintf('\n[Script 2 complete] Moody diagram with operating points plotted.\n');
