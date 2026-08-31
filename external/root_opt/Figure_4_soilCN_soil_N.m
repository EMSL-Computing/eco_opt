
clear all;
clc
close all
%% varying soil properties: Vmax_SOM and init_inorgN
terminal_time=120;
df_soilCN_inorgN_Vmax_SOM= struct();
init_inorgN = [10, 20];
CNSOM = [5,10,15,20,25,30]';CNSOM = (5:2.5:30)';
Vmax_SOM =[0.15,0.3];
numtstep=50;

root_col=[0.8350 0.3680 0.0000];  
exuflux_col = [0.0000 0.500 0.9]; 
supply_col = [1 1 1]*0.3;
leaf_col = [1 1 1]*0.5;
%%
[params, state_init]=params_base_condition();
for k =1:length(Vmax_SOM)
    for j= 1:length(init_inorgN)
        for i=1:length(CNSOM)
            params.Vmax_SOM= Vmax_SOM(k);
            state_init.inorgN=init_inorgN(j);
            state_init.SOMN = state_init.SOMC / CNSOM(i); % gN/m2
            [df, sol,~,~] = solve_ocp_nested(terminal_time, state_init, params,numtstep, []);
            fieldName = sprintf('df_VmaxSOM_%1.2f_Ni_%1.2f_SOMCN_%1.0f', Vmax_SOM(k),init_inorgN(j),CNSOM(i));
            fieldName = strrep(fieldName, '.', '_');
            df_soilCN_inorgN_Vmax_SOM.(fieldName) = df;
%             Lines  = readlines('out.txt');
%             exit_msg = Lines(32)
%             fig=plotting(df,params,exit_msg);
%             exportgraphics(fig, "figs/soil_property/"+fieldName+".png", Resolution=300)
%             close all
        end
    end
end
toc
save('df_soilCN_inorgN_Vmax_SOM.mat', 'df_soilCN_inorgN_Vmax_SOM')

%%
% soil soil figure with soil CN ratio on the X axis and variation with
% SOM decay rate and initial min N concentration

load('df_soilCN_inorgN_Vmax_SOM.mat')

title_font=13;
fig=figure;fig.Color='w';
fig.Position=[40,60,800,650];
[ax, ~] = tight_subplot(2,2, [.075 .1],[.1 .075],[0.1 .02]);

for w=1:4
    hold(ax(w),'on')
end
title(ax(1),"(A) Low min. N",'FontWeight','normal',FontSize=title_font);
title(ax(2),"(B) Low min. N",'FontWeight','normal',FontSize=title_font);
title(ax(3),"(C) High min. N",'FontWeight','normal',FontSize=title_font);
title(ax(4),"(D) High min. N",'FontWeight','normal',FontSize=title_font);

text(ax(2),-0.15, 0.25+0.4, 'or \int{\itG_{R}} / \int{\itS} [-]', 'Units', 'normalized', ...
    'Color', root_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);
text(ax(2), -0.15, 0.25, '{\it\intE / \intS} [-]', 'Units', 'normalized', ...
    'Color', exuflux_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);

text(ax(4),-0.15, 0.25+0.4, 'or \int{\itG_{R}} / \int{\itS} [-]', 'Units', 'normalized', ...
    'Color', root_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);
text(ax(4), -0.15, 0.25, '{\it\intE / \intS} [-]', 'Units', 'normalized', ...
    'Color', exuflux_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);

ylabel(ax(1),'\int{\itS} / \int{\itA_{net}} [-]')
ylabel(ax(3),'\int{\itS} / \int{\itA_{net}} [-]')

axgrid= [ax(1),ax(2);ax(3),ax(4)];
final_tbl=table();
w=1;
numstr =["(A)","(B)","(C)","(D)"];
ls = ["--","-"];msw=1;lw=2.0;title_font=14;

for j= 1:length(init_inorgN)
    for k =1:length(Vmax_SOM)
        Cpartition=[];
        for i=1:length(CNSOM)
            fieldName = sprintf('df_VmaxSOM_%1.2f_Ni_%1.2f_SOMCN_%1.0f', Vmax_SOM(k),init_inorgN(j),CNSOM(i));

            fieldName = strrep(fieldName, '.', '_');
            df=df_soilCN_inorgN_Vmax_SOM.(fieldName);
            S_by_Anet = trapz(df.time,df.rootCSupply)/trapz(df.time,df.Anet);
            Gleaf_by_Anet = trapz(df.time,df.leaf_growth_rate)/trapz(df.time,df.Anet);
            Groot_by_Anet = trapz(df.time,df.root_growth_rate)/trapz(df.time,df.Anet);
            exu_by_Anet = trapz(df.time,df.root_exu)/trapz(df.time,df.Anet);
            Groot_by_S = trapz(df.time,df.root_growth_rate)/trapz(df.time,df.rootCSupply);
            exu_by_S = trapz(df.time,df.root_exu)/trapz(df.time,df.rootCSupply);
            Cpartition = [Cpartition;[S_by_Anet,Gleaf_by_Anet,Groot_by_Anet, exu_by_Anet, Groot_by_S,   exu_by_S]];
        end
        tbl = array2table(Cpartition, "VariableNames",["S_by_Anet","Gleaf_by_Anet","Groot_by_Anet","exu_by_Anet","Gr_S","E_S"]);
        tbl.CNSOM=CNSOM;
        tbl.init_inorgN = repmat(init_inorgN(j), height(tbl), 1);
        tbl.Vmax_SOM = repmat(Vmax_SOM(k), height(tbl), 1);

        final_tbl=[final_tbl;tbl];

        plot(axgrid(j,1),tbl.CNSOM, tbl.S_by_Anet,LineStyle=ls(k), ...
            MarkerSize=msw, LineWidth=lw, Color=supply_col);
%         plot(axgrid(j,1),tbl.CNSOM, tbl.Gleaf_by_Anet,LineStyle=ls(k), ...
%             MarkerSize=msw, LineWidth=lw, Color=leaf_col);
        plot(axgrid(j,2),tbl.CNSOM, tbl.Gr_S,LineStyle=ls(k), ...
            MarkerSize=msw, LineWidth=lw, Color=root_col);
        plot(axgrid(j,2),tbl.CNSOM, tbl.E_S,LineStyle=ls(k), ...
            MarkerSize=msw, LineWidth=lw, Color=exuflux_col);
    end
end



% axis styling----------------------------

ax2 = axes;
ax2.Position = ax(1).Position;
p = [1,2];
for jj = 1:2
    p(jj) = plot(ax2, nan, nan, ls{jj}, 'LineWidth', 2, 'Color','k');
    hold(ax2, 'on')
end
hl = legend(ax2, p,  {'slow SOM decomposition','fast SOM decomposition'});
hl.FontSize = 11;hl.Location="southwest";
hl.Box = 'off';h1.Title.String="line color";
hl.Color='w';hl.LineWidth=0.5;
ax2.Visible = 'off';


% set(ax, 'Ylim',[0,0.4])
set(ax([1,3]), 'Ylim',[0,0.9])
set(ax([2,4]), 'Ylim',[0,0.4])

set(ax(1:2), 'Xticklabel',[])
% set(ax([2,4]), 'Yticklabel',[])

xlabel(ax(3:4),"SOM C:N ratio");

for i =1:length(ax)
    set(ax(i), 'LineWidth', 0.5, 'FontSize', 11, 'Box','on')
    ax(i).YLabel.FontSize=13;
    ax(i).XLabel.FontSize=13;
    ax(i).Title.FontSize=13;
    grid(ax(i),'on');
end

% exportgraphics(gcf, "figs/Figure3.png", Resolution=600)
% print(gcf, 'figs/Figure3.svg', '-dsvg');


%%
title_font = 13;
fig = figure; 
fig.Color = 'w';
fig.Position = [40, 60, 800, 800];

% Variables
varNames = ["Anet","S","Groot","exu"];

% Line styles for Vmax_SOM
ls = ["--","-"];  

% Create a 4x2 grid of axes
[ax, ~] = tight_subplot(4, 2, [.07 .07], [.07 .05], [.1 .04]);
for ii = 1:8
    hold(ax(ii), 'on');
end

final_tbl = table();

% -----------------------------
% MAIN LOOPS
% -----------------------------
for j = 1:length(init_inorgN)     % → COLUMN
    for k = 1:length(Vmax_SOM)    % → LINE STYLE

        Cpartition = nan(length(CNSOM), 4);

        for i = 1:length(CNSOM)
            fieldName = sprintf('df_VmaxSOM_%1.2f_Ni_%1.2f_SOMCN_%1.0f', ...
                                 Vmax_SOM(k), init_inorgN(j), CNSOM(i));
            fieldName = strrep(fieldName, '.', '_');
            df = df_soilCN_inorgN_Vmax_SOM.(fieldName);

            Cpartition(i, :) = [
                trapz(df.time, df.Anet), ...
                trapz(df.time, df.rootCSupply), ...
                trapz(df.time, df.root_growth_rate), ...
                trapz(df.time, df.root_exu)
            ];
        end

        tbl = array2table(Cpartition, "VariableNames", varNames);
        tbl.CNSOM      = CNSOM(:);
        tbl.init_inorgN = repmat(init_inorgN(j), height(tbl), 1);
        tbl.Vmax_SOM    = repmat(Vmax_SOM(k), height(tbl), 1);

        final_tbl = [final_tbl; tbl];

        % Plot (row = variable, column = init_inorgN)
        for v = 1:4
            row = v;
            col = j;
            idx = (row-1)*2 + col;

            plot(ax(idx), tbl.CNSOM, tbl.(varNames(v)), ...
                "LineStyle", ls(k), ...
                "Marker","o","MarkerSize",3,"MarkerEdgeColor",'k',...
                "Color", [1 1 1 ]*0.4, ...
                "LineWidth", 2);
        end
    end
end

% -----------------------------
% LABELS
% -----------------------------

ylabel(ax(1), "\int{{\itA_{net}} [gC m^{-2}]}", 'FontSize', 13, Interpreter='tex');
ylabel(ax(3), "\int{{\itS} [gC m^{-2}]}", 'FontSize', 13, Interpreter='tex');
ylabel(ax(5), "\int{{\itG_{root}} [gC m^{-2}]}", 'FontSize', 13, Interpreter='tex');
ylabel(ax(7), "\int{{\itE} [gC m^{-2}]}", 'FontSize', 13, Interpreter='tex');

ylim(ax([1,2]),[500,1600])
ylim(ax([3,4]),[400,850])
ylim(ax([5,6]),[100,300])
ylim(ax([7,8]),[0,120])

xlabel(ax(7), "SOM C:N ratio", 'FontSize', 13);
xlabel(ax(8), "SOM C:N ratio", 'FontSize', 13);

% -----------------------------
% TITLES FOR COLUMNS
% -----------------------------
title(ax(1),"Low min. N",'FontWeight','normal',FontSize=title_font);
title(ax(2),"High min. N",'FontWeight','normal',FontSize=title_font);
% -----------------------------
% LEGEND: Vmax line styles
% -----------------------------
ax_leg = axes;
ax_leg.Position = ax(1).Position;
hold(ax_leg, 'on');

p = gobjects(length(Vmax_SOM),1);
for k = 1:length(Vmax_SOM)
    p(k) = plot(ax_leg, nan, nan, ls(k), 'Color','k','LineWidth',2);
end

legend(ax_leg, p, {'slow SOM decomposition','fast SOM decomposition'}, ...
       "Location","northeast", "FontSize",11, "Box","off");

ax_leg.Visible = 'off';

% -----------------------------
% STYLE
% -----------------------------
for i = 1:length(ax)
    set(ax(i), 'LineWidth', 0.5, 'FontSize', 11, 'Box', 'on');
    grid(ax(i), 'on');
end

panel_labels = {'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)'};

for ii = 1:8
    text(ax(ii), 0.0, 1.1, panel_labels{ii}, ...
        'Units', 'normalized', ...
        'FontSize', 11);
end

% exportgraphics(fig, "figs/Figure3_SI.png", "Resolution", 600);



%%

close all
title_font = 16;
fig = figure; fig.Color='w';
fig.Position=[40,60,1400,700];

%-----------------------------------------------------------------------
% Create 2×4 grid
%-----------------------------------------------------------------------
[ax, ~] = tight_subplot(2,4,[.125 .06],[.1 .075],[.07 .03]);

for w=1:8
    hold(ax(w),'on')
end

for i=1:length(ax)
    set(ax(i),'LineWidth',0.5,'FontSize',title_font-3,'Box','on')
    grid(ax(i),'on')
    axis(ax(i), 'tight')    
end

% Panel titles (row 1 = low N, row 2 = high N)
title(ax(1), "(A) Low min. N",  'FontSize', title_font)
title(ax(2), "(B) Low min. N",  'FontSize', title_font)
title(ax(3), ["(C) Low min. N and", "{\it slow} SOM decomposition"],  'FontSize', title_font)
title(ax(4), ["(D) Low min. N and", "{\it fast} SOM decomposition"],  'FontSize', title_font)

title(ax(5), "(E) High min. N", 'FontSize', title_font)
title(ax(6), "(F) High min. N", 'FontSize', title_font)
title(ax(7), ["(G) High min. N and" ,"{\it slow} SOM decomposition"], 'FontSize', title_font)
title(ax(8), ["(H) High min. N and" ,"{\it fast} SOM decomposition"], 'FontSize', title_font)

ylabel(ax(1), '\int{\itS} / \int{\itA_{net}} [-]')
ylabel(ax(5), '\int{\itS} / \int{\itA_{net}} [-]')
xlabel(ax(5:8), "SOM C:N ratio")

ylabel(ax(1),'\int{\itS} / \int{\itA_{net}} [-]')
ylabel(ax(5),'\int{\itS} / \int{\itA_{net}} [-]')
xlabel(ax(5:8),"SOM C:N ratio")


text(ax(2),-0.2, 0.3+0.4, 'or \int{\itG_{R}} / \int{\itS} [-]', 'Units', 'normalized', ...
    'Color', root_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);
text(ax(2),-0.2, 0.2, '{\it\intE / \intS} [-]', 'Units', 'normalized', ...
    'Color', exuflux_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);

text(ax(6),-0.2, 0.3+0.4, 'or \int{\itG_{R}} / \int{\itS} [-]', 'Units', 'normalized', ...
    'Color', root_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);
text(ax(6),-0.2, 0.2, '{\it\intE / \intS} [-]', 'Units', 'normalized', ...
    'Color', exuflux_col, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', Rotation=90, FontSize=title_font);


ax_low  = ax(1:4);
ax_high = ax(5:8);

ls = ["--","-"]; lw=2;



%=======================================================================
% Storage to hold final leaf/root/exu vs CNSOM for stacked-area plots
%=======================================================================
leaf_final_low  = zeros(length(Vmax_SOM), length(CNSOM));
root_final_low  = zeros(length(Vmax_SOM), length(CNSOM));
exu_final_low   = zeros(length(Vmax_SOM), length(CNSOM));

leaf_final_high = zeros(length(Vmax_SOM), length(CNSOM));
root_final_high = zeros(length(Vmax_SOM), length(CNSOM));
exu_final_high  = zeros(length(Vmax_SOM), length(CNSOM));



%=======================================================================
% MAIN LOOP
%=======================================================================
for j = 1:length(init_inorgN)      % j=1 low N, j=2 high N

    for k = 1:length(Vmax_SOM)    % two Vmax

        Svals = []; GrSvals=[]; ESvals=[];
        leafF = []; rootF=[]; exuF=[];

        for i = 1:length(CNSOM)

            fieldName = sprintf('df_VmaxSOM_%1.2f_Ni_%1.2f_SOMCN_%1.0f', ...
                        Vmax_SOM(k), init_inorgN(j), CNSOM(i));
            fieldName = strrep(fieldName, '.', '_');
            df = df_soilCN_inorgN_Vmax_SOM.(fieldName);

            Svals(end+1)   = trapz(df.time,df.rootCSupply)      / trapz(df.time,df.Anet);
            GrSvals(end+1) = trapz(df.time,df.root_growth_rate) / trapz(df.time,df.rootCSupply);
            ESvals(end+1)  = trapz(df.time,df.root_exu)         / trapz(df.time,df.rootCSupply);

            leafF(end+1) = trapz(df.time,df.leaf_growth_rate)   ;
            rootF(end+1) = trapz(df.time,df.root_growth_rate);
            exuF(end+1)  = trapz(df.time,df.root_exu);
        end

        % Select correct row of axes
        if j == 1
            axrow = ax_low;
            leaf_final_low(k,:) = leafF;
            root_final_low(k,:) = rootF;
            exu_final_low(k,:)  = exuF;
        else
            axrow = ax_high;
            leaf_final_high(k,:) = leafF;
            root_final_high(k,:) = rootF;
            exu_final_high(k,:)  = exuF;
        end

        % ORIGINAL PLOTS (columns 1–2)
        plot(axrow(1), CNSOM, Svals,   LineStyle=ls(k), LineWidth=lw, Color=supply_col);
        plot(axrow(2), CNSOM, GrSvals, LineStyle=ls(k), LineWidth=lw, Color=root_col);
        plot(axrow(2), CNSOM, ESvals,  LineStyle=ls(k), LineWidth=lw, Color=exuflux_col);
    end
end



%=======================================================================
% STACKED AREA PLOTS vs CNSOM (col 3=Vmax1, col 4=Vmax2)
%=======================================================================
for k = 1:length(Vmax_SOM)

    % Row 1: low N
    area(ax_low(2+k), CNSOM, ...
        [leaf_final_low(k,:); root_final_low(k,:); exu_final_low(k,:)]')
    ax_low(2+k).Children(3).FaceColor = leaf_col;
    ax_low(2+k).Children(2).FaceColor = root_col;
    ax_low(2+k).Children(1).FaceColor = exuflux_col;
    ax_low(2+k).YLabel.String = '[gC m^{-2}]';

    % Row 2: high N
    area(ax_high(2+k), CNSOM, ...
        [leaf_final_high(k,:); root_final_high(k,:); exu_final_high(k,:)]')
    ax_high(2+k).Children(3).FaceColor = leaf_col;
    ax_high(2+k).Children(2).FaceColor = root_col;
    ax_high(2+k).Children(1).FaceColor = exuflux_col;
    ax_high(2+k).YLabel.String = '[gC m^{-2}]';
end

% -----------------------------
% LEGEND: Vmax line styles
% -----------------------------
ax_leg = axes;
ax_leg.Position = ax(1).Position;
hold(ax_leg, 'on');

p = gobjects(length(Vmax_SOM),1);
for k = 1:length(Vmax_SOM)
    p(k) = plot(ax_leg, nan, nan, ls(k), 'Color','k','LineWidth',2);
end

legend(ax_leg, p, {'slow SOM decomposition','fast SOM decomposition'}, ...
       "Location","southwest", "FontSize",13, "Box","off");
ax_leg.Visible = 'off';

% -----------------------------
% LEGEND: Area plot
% -----------------------------
legend(ax_low(3),{'\int{\it G_{L}}','\int{\it G_{R}}','\int{\it E}'}, 'Location','northeast', ...
    NumColumns=1, FontSize=14, Box='on')


%-----------------------------------------------------------------------
% Styling
%-----------------------------------------------------------------------

set(ax(1:4), 'Xticklabel',[])


set(ax([1,5]), 'Ylim',[0,0.9])
set(ax([2,6]), 'Ylim',[0,0.5]) 

set(ax([3,4,7,8]), 'Ylim',[0,1400]) 

exportgraphics(gcf,"figs/Figure4_revised.png",Resolution=600)
% print(gcf,'figs/Figure3.svg','-dsvg');



