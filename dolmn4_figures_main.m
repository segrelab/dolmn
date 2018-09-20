%REQUIRES CBREWER
% https://www.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2

clear variables
close all
clc

%% Load Data

if ~isdir('DOLMN_Plots/'); mkdir('DOLMN_Plots/'); end

loadDataName = 'Ecoli_iJR904';

% Plot Variables
legendSize = 16;
axesLabelSize = 16;
xyLabelSize = 20;
titleSize = 22;
lineWidth = 3;
markerSize = 100;
tickLength = 0.005;

% Axes Limits and Positions
trsptLim = [0,45];
intlLim = [200,285];
axPos = [0.1 0.1 0.82 0.82]; % axes position
cbarPos = [0.8 0.1 0.02 0.82]; % colorbar position

%% Figure 3

% Load Data
vars = {'trsptCon','intlCon','biomass1','biomass2','biomass1_nan','biomass2_nan','x1','y1','x2','y2'};
load(fullfile('DOLMN_Parsed',[loadDataName '_interp.mat']),vars{:})

% Colormaps
cmap_bioFlux = cbrewer('seq','Greens',256);
cmap_bioFluxDiff = cbrewer('seq','YlOrRd',256);

% Calculate Difference between 1- and 2-Strain Subnetworks
biomass_diff21 = biomass2 - biomass1;
biomass_diff21(biomass1 == 0) = NaN;
biomass_diff21(biomass_diff21 <= 0) = NaN;

% 3a: Biomass Flux for 1-Strain Subnetworks
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '3a'; ax = axes(fig);
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,biomass1_nan'); g.AlphaData = ~isnan(biomass1_nan'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFlux); caxis(ax,[0, 0.05*max(ceil(20.*[biomass1(:); biomass2(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 3a', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_3a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 3b: Biomass Flux for 2-Strain Subnetworks
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '3b'; ax = axes(fig);
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,biomass2_nan'); g.AlphaData = ~isnan(biomass2_nan'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFlux); caxis(ax,[0, 0.05*max(ceil(20.*[biomass1(:); biomass2(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 3b', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_3b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 3c: Difference in Biomass Flux for 1- and 2-Strain Subnetworks
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '3c'; ax = axes(fig);
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,biomass_diff21'); g.AlphaData = ~isnan(biomass_diff21'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFluxDiff); caxis(ax,[0, 0.05*max(ceil(20.*biomass_diff21(:)))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Difference in Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 3c', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_3c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure 4

% Load Data
vars1 = {'trsptCon','intlCon','jaccDist_2m','pathwayEucDist','pathwayJaccDist','pathwayNames','metNames','exchMets_2m','numExchMets_2m','x1','y1','x2','y2'};
load(fullfile('DOLMN_Parsed',[loadDataName '_interp.mat']),vars1{:})
vars2 = {'trsptCon2','score2m1','score2m2','explained','Nt2m','Ni2m'}; % PCA
load(fullfile('DOLMN_Parsed',[loadDataName '_summary.mat']),vars2{:})

% PCA Variables
Nt = [10, 19, 30, 45];
arrowHead = 10;
markerSize = 14;

% Colormaps
cmap_jaccDist = cbrewer('seq','Blues',256);
cmap_eucDist = cbrewer('seq','Reds',256);
cmap_exchMets = cbrewer('seq','Greys',nanmax(numExchMets_2m(:))+1);
cmap_succ = 0.45.*ones(1,3);
cmap_pca = flip(cbrewer('div','Spectral',numel(Nt)));

% 4a: Jaccard Distance for 2-Strain Subnetworks
n = 4;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '4a'; ax = axes(fig); clear h
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,jaccDist_2m'); g.AlphaData = ~isnan(jaccDist_2m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_jaccDist); caxis(ax,[0, nanmax(jaccDist_2m(:))]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 4a', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_4a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 4a (Inset): PCA of 2-Strain Subnetworks
n = 5;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '4a (Inset)'; ax = axes(fig); hold(ax,'on'); clear h
for ii = 1:numel(Nt)
    idx_t = find(Nt2m == Nt(ii));
    idx_i = find(Ni2m(idx_t) == 285);
    
    % Strain A
    x = score2m1(idx_t([idx_i,1]),1);
    y = score2m1(idx_t([idx_i,1]),2);
    p1 = [x(1),y(1)]; % first point
    mid_p = p1 + [diff(x),diff(y)]./2; % midway point
    v = [diff(x), diff(y)]; v(abs(v)<1E-6) = 0; % vector
    L = sqrt(sum(v.^2)); % vector length
    h(ii) = plot(ax,x,y,'-', 'LineWidth',1.5*lineWidth, 'Color',cmap_pca(ii,:), 'MarkerEdgeColor','none', 'MarkerFaceColor',cmap_pca(ii,:), 'MarkerSize',markerSize);
    ah = annotation(fig,'arrow', 'HeadStyle','vback1', 'HeadLength',arrowHead, 'HeadWidth',3*arrowHead, 'LineStyle','none');
    set(ah, 'Parent',gca, 'Position',[p1(1),p1(2),v(1),v(2)], 'Color',cmap_pca(ii,:))
    
    % Strain B
    x = score2m2(idx_t([idx_i,1]),1);
    y = score2m2(idx_t([idx_i,1]),2);
    p1 = [x(1),y(1)]; % first point
    mid_p = p1 + [diff(x),diff(y)]./2; % midway point
    v = [diff(x), diff(y)]; v(abs(v)<1E-6) = 0; % vector
    L = sqrt(sum(v.^2)); % vector length
    plot(ax,x,y,'--', 'LineWidth',1.5*lineWidth, 'Color',cmap_pca(ii,:), 'MarkerEdgeColor','none', 'MarkerFaceColor',cmap_pca(ii,:), 'MarkerSize',markerSize);
    ah = annotation(fig,'arrow', 'HeadStyle','vback1', 'HeadLength',arrowHead, 'HeadWidth',3*arrowHead, 'LineStyle','none');
    set(ah, 'Parent',gca, 'Position',[p1(1),p1(2),v(1),v(2)], 'Color',cmap_pca(ii,:))
end
hold(ax,'off'); box(ax,'on'); grid(ax,'on'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
ax.Position = axPos;
l = legend(h,cellstr(int2str(Nt'))); l.Location = 'SouthEast'; l.FontSize = legendSize; l.Box = 'off'; l.Title.String = 'T_{TR}';
xlabel(ax,['PC1: ' num2str(round(explained(1),1)) '%'], 'FontSize',xyLabelSize)
ylabel(ax,['PC2: ' num2str(round(explained(2),1)) '%'], 'FontSize',xyLabelSize)
title(ax,'Figure 4a (Inset)', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_4a_inset','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 4b: Number of Exchanged Metabolites for 2-Strain Subnetworks
n = 6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '4b'; ax = axes(fig); clear h
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,numExchMets_2m'); g.AlphaData = ~isnan(numExchMets_2m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_exchMets); caxis(ax,[0, nanmax(numExchMets_2m(:))+1]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (0:nanmax(numExchMets_2m(:)))+0.5; cbar.TickLabels = 0:nanmax(numExchMets_2m(:));
cbar.Label.String = 'Number of Exchanged Metabolites'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 4b', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_4b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 4c: Euclidean Distance for 2-Strain Subnetworks - TCA Cycle
[~,pathway_idx,~] = intersect(pathwayNames,'Citric Acid Cycle');
n = 7;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '4c'; ax = axes(fig); clear h
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,pathwayEucDist{pathway_idx}'); g.AlphaData = ~isnan(pathwayEucDist{pathway_idx}'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_eucDist); caxis(ax,[0, nanmax(pathwayEucDist{pathway_idx}(:))]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Euclidean Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 4c', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_4c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 4d: Where Succinate is Exchanged for 2-Strain Subnetworks
n = 8;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '4d'; ax = axes(fig); clear h
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
[~,met_idx,~] = intersect(metNames,'Succinate');
g = imagesc(ax,trsptCon,intlCon,exchMets_2m{met_idx}'); g.AlphaData = ~isnan(exchMets_2m{met_idx}'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_succ); caxis(ax,[0, 1]);
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure 4d', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_4d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars1{:})
clear(vars2{:})
clear cmap*

%% Figure 5

% Load Data
vars = {'rho2','P2_rho','exchMetsNum_2m','exchMetNames','mediumMets','metNames','biomass2'};
load(fullfile('DOLMN_Parsed',[loadDataName '_summary.mat']),vars{:})

% Colormaps
cmap_barPlot = 0.45.*ones(1,3);
cmap_rho = flip(cbrewer('div','RdBu',256));

% Hierarchical Clustering
rho = rho2(P2_rho,P2_rho);
exchMetNames = exchMetNames(P2_rho);
rm_idx = find(ismember(rho,zeros(1,size(rho,1)),'rows')); % remove all zeros
rho(rm_idx,:) = []; rho(:,rm_idx) = []; exchMetNames(rm_idx) = [];
[~,rm_idx,~] = intersect(exchMetNames,mediumMets); % remove medium metabolites
rho(rm_idx,:) = []; rho(:,rm_idx) = []; exchMetNames(rm_idx) = [];

% 5a: Bar Plot
n = 9;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '5a'; ax = axes(fig);
[~,idx_exchMet,idx_met] = intersect(exchMetNames,metNames,'stable');
N = 100.*exchMetsNum_2m(idx_met)./numel(find(biomass2 ~= 0));
h = bar(ax,N); h.BarWidth = 1; h.FaceColor = cmap_barPlot;
box(ax,'on'); grid(ax,'on'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
xlim(ax,[0.5,numel(idx_exchMet)+0.5])
ax.XTickLabelRotation = 90; ax.TickLabelInterpreter = 'none';
ax.XTick = 1:numel(exchMetNames); ax.XTickLabel = exchMetNames;
ax.Position = axPos;
ylabel(ax,'Percentage of Simulations Exchanged', 'FontSize',xyLabelSize)
title(ax,'Figure 5a', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_5a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% 5b: Spearman (Exchange Flux) Correlations
n = 10;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = '5b'; ax = axes(fig);
imagesc(ax,rho);
box(ax,'on'); grid(ax,'off'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
cbar = colorbar(ax); caxis(ax,[-1,1]); colormap(ax,cmap_rho);
cbar.Label.String = 'Spearman Correlation'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
ax.XTickLabelRotation = 90; ax.TickLabelInterpreter = 'none';
ax.XTick = 1:numel(exchMetNames); ax.XTickLabel = exchMetNames;
ax.YTick = 1:numel(exchMetNames); ax.YTickLabel = exchMetNames;
title(ax,'Figure 5b', 'FontSize',titleSize)
print('DOLMN_Plots/Figure_5b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars1{:})
clear(vars2{:})
clear cmap*
