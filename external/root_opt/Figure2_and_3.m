%% variation of initial soil inorg N and init CN ratio
close all
clear all
clc
%%

terminal_time=120;

df_t=table();
scenario=["base","low min. N","high SOM C:N"]; % if change the order here then change the df.scenario as well

% base scenario------------------------------------------------
[params, state_init]=params_base_condition();
[df, ~,~] = solve_ocp_nested(terminal_time, state_init, params,100, []);
Lines  = readlines('out.txt');
exit_msg = Lines(32);

% plotting(df,params,exit_msg);


df.scenario = repmat(scenario(1), height(df), 1);
df_t=[df_t;df];

% N_i scenario------------------------------------------------
[params, state_init]=params_base_condition();
state_init.inorgN = 7.5;
[df, ~,~] = solve_ocp_nested(terminal_time, state_init, params,100, []);

Lines  = readlines('out.txt');
exit_msg = Lines(32);
% plotting(df,params,exit_msg);
df.scenario = repmat(scenario(2), height(df), 1);
df_t=[df_t;df];

% soil CN scenario ------------------------------------------------
[params, state_init]=params_base_condition();
state_init.SOMN = state_init.SOMC /20; % gN/m2
[df, ~,~] = solve_ocp_nested(terminal_time, state_init, params,100, []);
Lines  = readlines('out.txt');
exit_msg = Lines(32);
% plotting(df,params,exit_msg);
df.scenario = repmat(scenario(3), height(df), 1);
df_t=[df_t;df];


%% Figure2_revised
close all

fig = figure;
fig.Position = [200 80 1100 1100];
fig.Color = 'w';

% 5 rows × 3 columns layout
t = tiledlayout(5,3);
t.TileSpacing = 'loose';
t.Padding = 'loose';
alpha_val = 0.95;
% -------------------------------------------------------------------------
% Create axes in required order
% -------------------------------------------------------------------------
axA = nexttile([1 2]);   % (A) A_net time series
axB = nexttile();        % (B) cumulative A_net

axC = nexttile([1 2]);   % (C) GL / A_net
axD = nexttile();        % (D) cumulative GL

axE = nexttile([1 2]);   % (E) S / A_net
axF = nexttile();        % (F) cumulative S

axG = nexttile([1 2]);   % (G) GR / S
axH = nexttile();        % (H) cumulative GR

axI = nexttile([1 2]);   % (I) E / S
axJ = nexttile();        % (J) cumulative E

ax = [axA axB axC axD axE axF axG axH axI axJ];

for k = 1:numel(ax)
    hold(ax(k),'on');
    box(ax(k),'on');
    grid(ax(k),'on');
end

color = copper(3);
% color = [
%     0.25 0.2500 0.25000   % grey
%     0.8350 0.3680 0.0000   % vermillion
%     0.0000 0.700 0.4510   % green
% ];
% -------------------------------------------------------------------------
% ---------- ROW 1: A_net time series and cumulative A_net -----------------
% -------------------------------------------------------------------------

for i = 1:length(scenario)
    df = df_t(strcmp(df_t.scenario, scenario(i)), :);
    cumA(i) = trapz(df.time, df.Anet);
    fprintf('%s  %.4f\n\n', scenario(i), cumA(i));
end

bar(axB, categorical(scenario), cumA, 'FaceAlpha',alpha_val, ...
    'FaceColor','flat','CData',color,'EdgeColor','flat');
scenario = string(axB.XTickLabel);

for i = 1:length(scenario)
    df = df_t(strcmp(df_t.scenario, scenario(i)), :);
    base_color = color(i,:);          % original RGB triplet
    line_color = base_color * alpha_val + (1 - alpha_val) * [1 1 1];  
    % (A) time series
    plot(axA, df.time, df.Anet, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
end



ylabel(axA, '{\it A_{net}} [gC m^{-2} d^{-1}]');
ylabel(axB, '{\it \int A_{net}} [gC m^{-2}]');


% -------------------------------------------------------------------------
% ---------- COMMON VARIABLES FOR ROWS 2–5 --------------------------------
% -------------------------------------------------------------------------
yvar_ratio = ["leaf_growth_rate","rootCSupply","root_growth_rate","root_exu"];
ratio_den  = ["Anet","Anet","Net_rootC_assimilation","Net_rootC_assimilation"];
ylabels_ratio = ["G_{L}/A_{net}", "S/A_{net}", "G_{R}/S_{net}", "E/S_{net}"];
ylims_ratio = [1,1, 1, 0.6];  

ax_ratio = [axC axE axG axI];
ax_bar =   [axD axF axH axJ];
 
% -------------------------------------------------------------------------
% ---------- ROWS 2–5 -----------------------------------------------------
% -------------------------------------------------------------------------

for j = 1:4
    temp = [];
    for i = 1:length(scenario)
        df = df_t(strcmp(df_t.scenario, scenario(i)), :);
        ynum = df.(yvar_ratio(j));
        % Cumulative
        temp(i) = trapz(df.time, ynum);
        fprintf('%s %s  %.4f\n\n', scenario(i), yvar_ratio(j), temp(i));
    end
    % Cumulative bar chart
    b = bar(ax_bar(j), categorical(scenario), temp);
    b.FaceColor = 'flat';
    b.CData = color;
    b.FaceAlpha = alpha_val;
    b.EdgeColor = 'flat';
    ylabel(ax_bar(j), "\int " + extractBefore(ylabels_ratio(j),"/") + " [gC m^{-2}]");

end
scenario = string(ax_bar(j).XTickLabel);


for j = 1:4
    temp = [];
    for i = 1:length(scenario)
        df = df_t(strcmp(df_t.scenario, scenario(i)), :);
        % Ratio time series
        ynum = df.(yvar_ratio(j));
        yden = df.(ratio_den(j));
        base_color = color(i,:);          % original RGB triplet
        line_color = base_color * alpha_val + (1 - alpha_val) * [1 1 1];  
        plot(ax_ratio(j), df.time, ynum ./ yden, ...
            'LineWidth', 2, 'Color', line_color, 'DisplayName', scenario(i));

        % Cumulative
        temp(i) = trapz(df.time, ynum);
    end
    ylabel(ax_ratio(j), "{\it" + ylabels_ratio(j) + "} [-]");
    ylim(ax_ratio(j), [0, ylims_ratio(j)]);
end

yline(axC, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axE, 0.8, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axG, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axG, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)

yline(axI, params.min_exudation , '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axI, params.max_exudation , '--', 'Color',[1 1 1]*0.25, LineWidth=1)

set(axE, 'Ylim',[0, 1])
set([axD,axF], 'Ylim',[0, 900])
set([axH, axJ], 'Ylim',[0, 300])

% -------------------------------------------------------------------------
% Formatting
% -------------------------------------------------------------------------
% Remove xticklabels on all except bottom row
set([axA,axB], 'XTickLabel', []);
set([axC,axD], 'XTickLabel', []);
set([axE,axF], 'XTickLabel', []);
set([axG,axH], 'XTickLabel', []);

xlabel(axI, 'time [d]');

lh=legend(axA);
lh.Location="northwest";
lh.FontSize=14;lh.NumColumns=2;lh.Box='off';
% -------------------------------------------------------------------------
% Titles (A to J)
% -------------------------------------------------------------------------
titles = ["A","B","C","D","E","F","G","H","I","J"];
for k = 1:numel(ax)
    ttl = title(ax(k), "(" + titles(k) + ")", 'FontWeight','normal');
    ttl.Units = 'normalized';
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';
    ttl.FontSize = 16;

    set(ax(k), 'LineWidth', 0.5, 'FontSize',14, 'Box','on')
    ax(k).YLabel.FontSize=16;
    ax(k).XLabel.FontSize=16;
    grid(ax(k),'on');
end

exportgraphics(fig, "figs/Figure2_revised.png", Resolution=600);
% print(fig, 'figs/Figure2.svg', '-dsvg');


%% Figure2_revised_v2

fig = figure;
fig.Position = [200 80 1100 1100];
fig.Color = 'w';

% 5 rows × 3 columns layout
t = tiledlayout(5,3);
t.TileSpacing = 'loose';
t.Padding = 'loose';
alpha_val = 0.95;
% -------------------------------------------------------------------------
% Create axes in required order
% -------------------------------------------------------------------------
axA = nexttile([1 2]);   % (A) A_net time series
axB = nexttile();        % (B) cumulative A_net

axC = nexttile([1 2]);   % (C) GL / A_net
axD = nexttile();        % (D) cumulative GL

axE = nexttile([1 2]);   % (E) S / A_net
axF = nexttile();        % (F) cumulative S

axG = nexttile([1 2]);   % (G) GR / S
axH = nexttile();        % (H) cumulative GR

axI = nexttile([1 2]);   % (I) E / S
axJ = nexttile();        % (J) cumulative E

ax = [axA axB axC axD axE axF axG axH axI axJ];

for k = 1:numel(ax)
    hold(ax(k),'on');
    box(ax(k),'on');
    grid(ax(k),'on');
end

% color = copper(3);
% color = [
%     0.25 0.2500 0.25000   % grey
%     0.8350 0.3680 0.0000   % vermillion
%     0.0000 0.700 0.4510   % bluish-green
% ];

% -------------------------------------------------------------------------
% ---------- ROW 1: A_net time series and cumulative A_net -----------------
% -------------------------------------------------------------------------
for i = 1:length(scenario)
    df = df_t(strcmp(df_t.scenario, scenario(i)), :);
    base_color = color(i,:);          % original RGB triplet
    line_color = base_color * alpha_val + (1 - alpha_val) * [1 1 1];  
    % (A) time series
    plot(axA, df.time, df.Anet, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));

    % (B) cumulative
    cumA(i) = trapz(df.time, df.Anet);
    fprintf('%s  %.4f\n', scenario(i), cumA(i));
end

bar(axB, categorical(scenario), cumA, 'FaceAlpha',alpha_val, ...
    'FaceColor','flat','CData',color,'EdgeColor','flat');

ylabel(axA, '{\it A_{net}} [gC m^{-2} d^{-1}]');
ylabel(axB, '{\it \int A_{net}} [gC m^{-2}]');


% -------------------------------------------------------------------------
% ---------- COMMON VARIABLES FOR ROWS 2–5 --------------------------------
% -------------------------------------------------------------------------
yvar_ratio = ["leaf_growth_rate","rootCSupply","root_growth_rate","root_exu"];
ratio_den  = ["Anet","Anet","Net_rootC_assimilation","Net_rootC_assimilation"];
ylabels_ratio = ["G_{L}/A_{net}", "S/A_{net}", "G_{R}/S_{net}", "E/S_{net}"];

ylabels_ratio_bar = ["\int G_{L}/\int A_{net}", "\int S/\int A_{net}", "\int G_{R}/\int S_{net}", "\int E/\int S_{net}"];

ylims_ratio = [1,1, 1, 0.6];  

ax_ratio = [axC axE axG axI];
ax_bar =   [axD axF axH axJ];
 
% -------------------------------------------------------------------------
% ---------- ROWS 2–5 -----------------------------------------------------
% -------------------------------------------------------------------------
for j = 1:4

    temp = [];

    for i = 1:length(scenario)
        df = df_t(strcmp(df_t.scenario, scenario(i)), :);

        % Ratio time series
        ynum = df.(yvar_ratio(j));
        yden = df.(ratio_den(j));
        base_color = color(i,:);          % original RGB triplet
        line_color = base_color * alpha_val + (1 - alpha_val) * [1 1 1];  
        plot(ax_ratio(j), df.time, ynum ./ yden, ...
            'LineWidth', 2, 'Color', line_color, 'DisplayName', scenario(i));

        % Cumulative
        %temp(i) = trapz(df.time, ynum);
        temp(i) = trapz(df.time, ynum)/trapz(df.time, yden);
        fprintf('%s %s  %.4f\n', scenario(i), yvar_ratio(j), temp(i));

    end


    % Cumulative bar chart
    b = bar(ax_bar(j), categorical(scenario), temp);
    b.FaceColor = 'flat';
    b.CData = color;
    b.FaceAlpha = alpha_val;
    b.EdgeColor = 'flat';

    ylabel(ax_ratio(j), "{\it" + ylabels_ratio(j) + "} [-]");
%     ylabel(ax_bar(j), "\int " + extractBefore(ylabels_ratio(j),"/") + " [gC m^{-2}]");
    ylabel(ax_bar(j), ylabels_ratio_bar(j));


    ylim(ax_ratio(j), [0, ylims_ratio(j)]);
end

yline(axC, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axE, 0.8, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axG, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axG, 0.2, '--', 'Color',[1 1 1]*0.25, LineWidth=1)

yline(axI, params.min_exudation , '--', 'Color',[1 1 1]*0.25, LineWidth=1)
yline(axI, params.max_exudation , '--', 'Color',[1 1 1]*0.25, LineWidth=1)

set(axE, 'Ylim',[0, 1])
% Set limits and tick positions
set([axD, axF], 'YLim', [0 0.75], ...
                'YTick', [0 0.25 0.5 0.75]);
set([axH, axJ], 'YLim', [0 0.75], ...
                'YTick', [0 0.25 0.5 0.75]);

% set([axD,axF], 'Ylim',[0, 900])
% set([axH, axJ], 'Ylim',[0, 300])

% -------------------------------------------------------------------------
% Formatting
% -------------------------------------------------------------------------
% Remove xticklabels on all except bottom row
set([axA,axB], 'XTickLabel', []);
set([axC,axD], 'XTickLabel', []);
set([axE,axF], 'XTickLabel', []);
set([axG,axH], 'XTickLabel', []);

xlabel(axI, 'time [d]');

% -------------------------------------------------------------------------
% Titles (A to J)
% -------------------------------------------------------------------------
titles = ["A","B","C","D","E","F","G","H","I","J"];
for k = 1:numel(ax)
    ttl = title(ax(k), "(" + titles(k) + ")", 'FontWeight','normal');
    ttl.Units = 'normalized';
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';
    ttl.FontSize = 16;

    set(ax(k), 'LineWidth', 0.5, 'FontSize',14, 'Box','on')
    ax(k).YLabel.FontSize=16;
    ax(k).XLabel.FontSize=16;
    grid(ax(k),'on');
end

% exportgraphics(fig, "figs/Figure2_revised_v2.png", Resolution=600);


%% Figure2_SI 

fig = figure;fig.Position=[224   100   900   900];
fig.Color='w';
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
[ax, pos] = tight_subplot(5,3, [.02,.085],[.05 .015],[0.1 .02]);
for i=1:15
    hold(ax(i),'on');
end
lstyle = ["-","-","-"];
lw = [1,1,1]*2;
ylabel(ax(1),'leaf C [gC m^{-2}]')
ylabel(ax(2),'root C [gC m^{-2}]')
ylabel(ax(3),'exudate C [gC m^{-2}]')
ylabel(ax(4),'SOM C [gC m^{-2}]')
ylabel(ax(5),'mic C [gC m^{-2}]')
ylabel(ax(6),'min. N [gN m^{-2}]')
ylabel(ax(7),["net photosyn. rate","{\itA_{net}} [gC m^{-2} d^{-1}]"])
ylabel(ax(8),["root C supply rate","{\itS} [gC m^{-2} d^{-1}]"])
ylabel(ax(9),["leaf growth rate","{\itG_{L}} [gC m^{-2} d^{-1}]"])
ylabel(ax(10),["root growth rate","{\itG_{R}} [gC m^{-2} d^{-1}]"])
ylabel(ax(11),["root exu. rate","{\itE} [gC m^{-2} d^{-1}]"])
ylabel(ax(12),["tot. root resp. rate","[gC m^{-2} d^{-1}]"])
ylabel(ax(13),["root N uptake rate","{\itU_N} [gN m^{-2} d^{-1}]"])
ylabel(ax(14),["leaf N demand","{\it p_N} [gN m^{-2} d^{-1}]"])
ylabel(ax(15),"U_N/U_{N, pot} [-]")

color=copper(3);
for i =1:3
    df = df_t(strcmp(df_t.scenario,scenario(i)), :);
    plot(ax(1), df.time, df.leafC, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(2),df.time, df.rootC, 'linewidth', lw(i), ...
        'LineStyle' ,lstyle(i),"Color",color(i,:));

    plot(ax(3),df.time, df.exuC, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(4),df.time, df.SOMC, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(5),df.time, df.micC, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(6),df.time, df.inorgN, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(7),df.time, df.Anet, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));
    
    plot(ax(8),df.time, df.rootCSupply, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));
    
    plot(ax(9),df.time, df.leaf_growth_rate, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(10),df.time, df.root_growth_rate, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(11),df.time, df.root_exu, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(12),df.time, df.root_resp+df.root_maintenance_resp, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(13),df.time, df.root_N_uptake, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(14),df.time, df.leaf_N_uptake, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));

    plot(ax(15),df.time, df.rootN_uptake_norm, 'linewidth', lw(i), ...
        'LineStyle',lstyle(i),  "Color",color(i,:));
end

% set(ax([1:4,6]), 'Ylim',[0, inf])
set(ax(15), 'Ylim',[0, 1])
set(ax(1:12), 'Xticklabel',[])
xlabel(ax(13:15),'time [d]');

lh=legend(ax(4), scenario);
lh.Title.String="{Model scenarios}";lh.Location="southwest";
lh.FontSize=9;lh.NumColumns=1;lh.Box='on';

strs='ABCDEFGHIJKLMNOPQRST';
for i =1:length(ax)
    set(ax(i), 'LineWidth', 0.5, 'FontSize',10.5, 'Box','on')
%     ax(i).YLabel.FontSize=12;
%     ax(i).XLabel.FontSize=12;
    grid(ax(i),'on');
%     ttl = title(ax(i),"("+strs(i)+")",'FontWeight','normal',FontSize=10);
%     ttl.Units = 'Normalize';
%     ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
%     ttl.HorizontalAlignment = 'center';
end


% exportgraphics(gcf, "figs/Figure2_SI.png", Resolution=300)
% print(gcf, 'figs/Figure2_SI.svg', '-dsvg');



%% Figure - N uptake dynamics and controls

fig = figure;
fig.Position = [100 100 1000 900];
fig.Color = 'w';

t = tiledlayout(4, 2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Create axes
ax1 = nexttile([1 2]);   % Row 1: rootN_uptake_norm (spans 2 columns)
ax2 = nexttile();        % Row 2 Col 1: root growth rate
ax3 = nexttile();        % Row 2 Col 2: root N uptake rate
ax4 = nexttile();        % Row 3 Col 1: exudation rate
ax5 = nexttile();        % Row 3 Col 2: phiN
ax6 = nexttile();        % Row 4 Col 1: UN vs inorgN
ax7 = nexttile();        % Row 4 Col 2: UN vs rootN

ax = [ax1 ax2 ax3 ax4 ax5 ax6 ax7];

for i = 1:numel(ax)
    hold(ax(i), 'on');
    box(ax(i), 'on');
    grid(ax(i), 'on');
end

% Define colors matching other figures
color = copper(3);

% Loop through scenarios
for i = 1:length(scenario)
    df = df_t(strcmp(df_t.scenario, scenario(i)), :);
    
    % Panel 1: rootN_uptake_norm vs time (top, spanning 2 columns)
    plot(ax1, df.time, df.rootN_uptake_norm, ...
        'LineWidth', 2.5, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 2: root growth rate vs time
    plot(ax2, df.time, df.root_growth_rate, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 3: exudation rate vs time
    plot(ax3, df.time, df.root_exu, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 4: root N uptake rate vs time
    plot(ax4, df.time, df.root_N_uptake, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 5: phiN vs time
    plot(ax5, df.time, df.phiN, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 6: UN vs inorgN
    plot(ax6, df.inorgN, df.root_N_uptake, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    % Panel 7: UN vs rootN
    plot(ax7, df.rootN, df.root_N_uptake, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
end

% Labels
xlabel(ax1, 'time [d]');
ylabel(ax1, 'U_N/U_{N,pot} [-]');

xlabel(ax2, 'time [d]');
ylabel(ax2, {'root growth rate', '{\itG_R} [gC m^{-2} d^{-1}]'});

xlabel(ax3, 'time [d]');
ylabel(ax3, {'exudation rate', '{\itE} [gC m^{-2} d^{-1}]'});

xlabel(ax4, 'time [d]');
ylabel(ax4, {'root N uptake rate', '{\itU_N} [gN m^{-2} d^{-1}]'});

xlabel(ax5, 'time [d]');
ylabel(ax5, {'net mineralisation', '\phi_N [gN m^{-2} d^{-1}]'});

xlabel(ax6, 'inorg N [gN m^{-2}]');
ylabel(ax6, '{\itU_N} [gN m^{-2} d^{-1}]');

xlabel(ax7, 'root N [gN m^{-2}]');
ylabel(ax7, '{\itU_N} [gN m^{-2} d^{-1}]');

% Set specific limits
ylim(ax1, [0, 1.1]);

% Remove x-tick labels for top row
set(ax1, 'XTickLabel', []);

% Add panel labels
panel_labels = ["A", "B", "C", "D", "E", "F", "G"];
for i = 1:numel(ax)
    ttl = title(ax(i), "(" + panel_labels(i) + ")", 'FontWeight', 'normal');
    ttl.Units = 'normalized';
    ttl.Position(1) = 0;
    ttl.HorizontalAlignment = 'left';
    ttl.FontSize = 14;
    
    set(ax(i), 'LineWidth', 0.5, 'FontSize', 11);
    ax(i).YLabel.FontSize = 13;
    ax(i).XLabel.FontSize = 13;
end

% Add legend to top panel
lh = legend(ax1, 'Location', 'best', 'Orientation', 'horizontal');
lh.FontSize = 11;
lh.Box = 'on';

% exportgraphics(fig, "figs/Figure2_SI_revised.png", Resolution=600);
% print(fig, 'figs/Figure2_N_uptake_controls.svg', '-dsvg');

%%



fig = figure;
fig.Position = [100 100 700 850];
fig.Color = 'w';

t = tiledlayout(4, 1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Create axes
ax2 = nexttile();        % Row 1 Col 1: phiN  
ax3 = nexttile();        % Row 1 Col 2: root N uptake rate
ax4 = nexttile();        % Row 2 Col 1: root growth rate 
ax5 = nexttile();        % Row 2 Col 2: exudation rate 


ax = [ax2 ax3 ax4 ax5];

for i = 1:numel(ax)
    hold(ax(i), 'on');
    box(ax(i), 'on');
    grid(ax(i), 'on');
end

% Define colors matching other figures
color = copper(3);

% Loop through scenarios
for i = 1:length(scenario)
    df = df_t(strcmp(df_t.scenario, scenario(i)), :);
    

    plot(ax2, df.time, df.root_N_uptake  , ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    plot(ax3, df.time, df.root_growth_rate, ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
    plot(ax4, df.time, df.root_exu , ...
        'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));

    plot(ax5, df.time, df.phiN, ...
    'LineWidth', 2, 'Color', color(i,:), 'DisplayName', scenario(i));
    
end


ylabel(ax3, {'root growth rate', '{\itG_R} [gC m^{-2} d^{-1}]'});

ylabel(ax4, {'exudation rate', '{\itE} [gC m^{-2} d^{-1}]'});

ylabel(ax2, {'root N uptake rate', '{\itU_N} [gN m^{-2} d^{-1}]'});

xlabel(ax5, 'time [d]');
ylabel(ax5, {'net mineralisation', '\phi_N [gN m^{-2} d^{-1}]'});

xticklabels([ax2,ax3,ax4],[])


% Add panel labels
panel_labels = ["A", "B", "C", "D"];
for i = 1:numel(ax)
    ttl = title(ax(i), "(" + panel_labels(i) + ")", 'FontWeight', 'normal');
    ttl.Units = 'normalized';
    ttl.Position(1) = 0;
    ttl.HorizontalAlignment = 'left';
    ttl.FontSize = 14;
    
    set(ax(i), 'LineWidth', 0.5, 'FontSize', 12);
    ax(i).YLabel.FontSize = 13;
    ax(i).XLabel.FontSize = 13;
end

% Add legend to top panel
lh = legend(ax2, 'Location', 'best', 'Orientation', 'horizontal');
lh.FontSize = 14;
lh.Box = 'on';

exportgraphics(fig, "figs/Figure3.png", Resolution=600);
% print(fig, 'figs/Figure2_N_uptake_controls.svg', '-dsvg');







