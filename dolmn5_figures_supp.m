%REQUIRES CBREWER
%https://www.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2

clear variables
close all
clc

%% Load Data

if ~isdir('DOLMN_Plots_Supp/'); mkdir('DOLMN_Plots_Supp/'); end

loadDataName_full = 'Ecoli_iJR904';
loadDataName_core = 'Ecoli_Core';

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

%% Figure S1: Minimum number of intracellular reactions needed for growth

% Load Data
vars1a = {'intlCon*','intlRxns2m_same','intlRxns3m_same3'}; vars1b = {'minNi'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars1a{:})
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars1b{:})

% Colormaps
cmap_2m = [230, 85, 13]./256;
cmap_3m = [117,107,177]./256;
cmap_rxns2m = cbrewer('seq','Oranges',256);
cmap_rxns3m = cbrewer('seq','Purples',256);

% S1a: Min Ni
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S1a'; ax = axes(fig); clear h*
plot(ax,1:3,minNi,'ks-', 'LineWidth',lineWidth, 'MarkerFaceColor','k')
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
box(ax,'on'); grid(ax,'on'); xlim(ax,[0,4]); ylim(ax,intlLim); ax.XTick = 1:3;
xlabel(ax,'K', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S1a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S1a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

[~,intlIdx_2m,~] = find(intlCon2 == minNi(2));
[~,intlIdx_3m,~] = find(intlCon3 == minNi(3));
% S1b: Same Ni
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S1b'; ax = axes(fig); clear h*
h(1) = errorbar(ax,2,nanmean(intlRxns2m_same(:,intlIdx_2m)),nanstd(intlRxns2m_same(:,intlIdx_2m)),'s', 'LineWidth',lineWidth); hold(ax,'on')
h(1).Color = cmap_2m; h(1).MarkerFaceColor = cmap_2m;
h(2) = errorbar(ax,3,nanmean(intlRxns3m_same3(:,intlIdx_3m)),nanstd(intlRxns3m_same3(:,intlIdx_3m)),'s', 'LineWidth',lineWidth); hold(ax,'off')
h(2).Color = cmap_3m; h(2).MarkerFaceColor = cmap_3m;
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
box(ax,'on'); grid(ax,'on'); xlim(ax,[0,4]); ylim(ax,[150,165]); ax.XTick = 1:3;
legend(h,{'2 Strains','3 Strains'}, 'Location','NorthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'K', 'FontSize',xyLabelSize)
ylabel(ax,'# of Intracellular Reactions the Same', 'FontSize',xyLabelSize)
title(ax,'Figure S1b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S1b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

clear(vars1a{:}); clear intlIdx_2m intlIdx_3m
vars2 = {'trsptCon','intlCon','intlRxns2m_same','intlRxns3m_same3','x1','y1','x2','y2','x3','y3'};
load(fullfile(pwd,'DOLMN_Parsed',[loadDataName_full '_interp.mat']),vars2{:})

% S1c: Unique Intracellular (2 Models)
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S1c'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,intlRxns2m_same'); g.AlphaData = ~isnan(intlRxns2m_same'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_rxns2m); caxis(ax,[0, max([nanmax(intlRxns2m_same(:)),nanmax(intlRxns3m_same3(:)),intlLim])])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Number of Unique Intracellular Reactions'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S1c', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S1c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S1d: Unique Intracellular (3 Models)
n = 4;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S1d'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,intlRxns3m_same3'); g.AlphaData = ~isnan(intlRxns3m_same3'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_rxns3m); caxis(ax,[0, max([nanmax(intlRxns2m_same(:)),nanmax(intlRxns3m_same3(:)),intlLim])])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Number of Unique Intracellular Reactions'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S1d', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S1d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars1b{:}); clear(vars2{:})
clear cmap* vars*

%% Figure S2: (Nt,Ni) landscapes for Core 2-strain subnetworks

% Load Data
vars = {'trsptCon','intlCon','biomass*','jaccDist_2m','numExchMets_2m','tcaEucDist','x1','y1','x2','y2'};
load(fullfile('DOLMN_Parsed',[loadDataName_core '_interp.mat']),vars{:})

% Colormaps
cmap_bioFlux = cbrewer('seq','Greens',256);
cmap_bioFluxDiff = cbrewer('seq','YlOrRd',256);
cmap_jaccDist = cbrewer('seq','Blues',256);
cmap_exchMets = cbrewer('seq','Greys',nanmax(numExchMets_2m(:)+1));

% Calculate Differences
% 1- and 2-Strain Subnetworks
biomass_diff21 = biomass2 - biomass1;
biomass_diff21(biomass1 == 0) = NaN;
biomass_diff21(biomass_diff21 <= 0) = NaN;

% S2a: Biomass Flux for 1-Strain Subnetworks
n = 5;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S2a'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,biomass1_nan'); g.AlphaData = ~isnan(biomass1_nan'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0,25, 0,50],'square');
colormap(ax,cmap_bioFlux); caxis(ax,[0, 0.05*max(ceil(20.*[biomass1(:); biomass2(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S2a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S2a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S2b: Biomass Flux for 2-Strain Subnetworks
n = 6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S2b'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,biomass2_nan'); g.AlphaData = ~isnan(biomass2_nan'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0,25, 0,50],'square');
colormap(ax,cmap_bioFlux); caxis(ax,[0, 0.05*max(ceil(20.*[biomass1(:); biomass2(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S2b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S2b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S2c: Jaccard Distance for 2-Strain Subnetworks
n = 7;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S2c'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,jaccDist_2m'); g.AlphaData = ~isnan(jaccDist_2m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0,25, 0,50],'square');
colormap(ax,cmap_jaccDist); caxis(ax,[0, 0.05*nanmax(ceil(20.*jaccDist_2m(:)))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S2c', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S2c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S2d: Number of Exchanged Metabolites for 3-Strain Subnetworks
n = 8;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S2d'; ax = axes(fig); clear h*
patch(ax,[0,6.5,6.5,0],[0,0,50,50],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,numExchMets_2m'); g.AlphaData = ~isnan(numExchMets_2m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0,25, 0,50],'square');
colormap(ax,cmap_exchMets); caxis(ax,[0, nanmax(numExchMets_2m(:))+1])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (0:nanmax(numExchMets_2m(:)))+0.5; cbar.TickLabels = 0:nanmax(numExchMets_2m(:));
cbar.Label.String = 'Number of Exchanged Metabolites'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S2d', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S2d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S3: (Nt,Ni) landscapes for 3-strain subnetworks

% Load Data
vars = {'trsptCon','intlCon','biomass*','jaccDist*','numExchMets*','x1','y1','x2','y2','x3','y3'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_interp.mat']),vars{:})

% Colormaps
cmap_bioFlux = cbrewer('seq','Greens',256);
cmap_bioFluxDiff = cbrewer('seq','YlOrRd',256);
cmap_jaccDist = cbrewer('seq','Blues',256);
cmap_exchMets = cbrewer('seq','Greys',nanmax(numExchMets_3m(:))+1);

% Calculate Differences
% 1- and 2-Strain Subnetworks
biomass_diff21 = biomass2 - biomass1;
biomass_diff21(biomass1 == 0) = NaN;
biomass_diff21(biomass_diff21 <= 0) = NaN;
% 1- and 3-Strain Subnetworks
biomass_diff31 = biomass3 - biomass1;
biomass_diff31(biomass1 == 0) = NaN;
biomass_diff31(biomass_diff31 <= 0) = NaN;
% 2- and 3-Strain Subnetworks
biomass_diff32 = biomass3 - biomass2;
biomass_diff32(biomass2 == 0) = NaN;
biomass_diff32(biomass_diff32 <= 0) = NaN;

% S3a: Biomass Flux for 3-Strain Subnetworks
n = 9;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3a'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,biomass3_nan'); g.AlphaData = ~isnan(biomass3_nan'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFlux); caxis(ax,[0, 0.05*max(ceil(20.*[biomass1(:); biomass2(:); biomass3(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3b: Difference in Biomass Flux for 1- and 3-Strain Subnetworks
n = 10;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3b'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,biomass_diff31'); g.AlphaData = ~isnan(biomass_diff31'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFluxDiff); caxis(ax,[0, 0.05*max(ceil(20.*[biomass_diff21(:); biomass_diff31(:); biomass_diff32(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Difference in Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
legend(h,{'1 Strain'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
ax.Position = axPos; cbar.Position = cbarPos;
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3c: Difference in Biomass Flux for 2- and 3-Strain Subnetworks
n = 11;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3c'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,biomass_diff32'); g.AlphaData = ~isnan(biomass_diff32'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[trsptLim, intlLim],'square');
colormap(ax,cmap_bioFluxDiff); caxis(ax,[0, 0.05*max(ceil(20.*[biomass_diff21(:); biomass_diff31(:); biomass_diff32(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Difference in Biomass Flux (hr^{-1})'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3c', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3d: Jaccard Distance for 3-Strain Subnetworks, Strains A & B
n = 12;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3d'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,jaccDist_3m12'); g.AlphaData = ~isnan(jaccDist_3m12'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_jaccDist); caxis(ax,[0, 0.05*max(ceil(20.*[jaccDist_3m(:); jaccDist_3m12(:); jaccDist_3m13(:); jaccDist_3m23(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3d', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3e: Jaccard Distance for 3-Strain Subnetworks, Strains A & C
n = 13;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3e'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,jaccDist_3m13'); g.AlphaData = ~isnan(jaccDist_3m13'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_jaccDist); caxis(ax,[0, 0.05*max(ceil(20.*[jaccDist_3m(:); jaccDist_3m12(:); jaccDist_3m13(:); jaccDist_3m23(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3e', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3e','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3f: Jaccard Distance for 3-Strain Subnetworks, Strains B & C
n = 14;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3f'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,jaccDist_3m23'); g.AlphaData = ~isnan(jaccDist_3m23'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_jaccDist); caxis(ax,[0, 0.05*max(ceil(20.*[jaccDist_3m(:); jaccDist_3m12(:); jaccDist_3m13(:); jaccDist_3m23(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3f', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3f','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3g: Jaccard Distance for 3-Strain Subnetworks, Average
n = 15;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3g'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,jaccDist_3m'); g.AlphaData = ~isnan(jaccDist_3m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_jaccDist); caxis(ax,[0, 0.05*max(ceil(20.*[jaccDist_3m(:); jaccDist_3m12(:); jaccDist_3m13(:); jaccDist_3m23(:)]))])
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3g', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3g','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S3h: Number of Exchanged Metabolites for 3-Strain Subnetworks
n = 16;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S3h'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
g = imagesc(ax,trsptCon,intlCon,numExchMets_3m'); g.AlphaData = ~isnan(numExchMets_3m'); clear g
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
colormap(ax,cmap_exchMets); caxis(ax,[0, nanmax(numExchMets_3m(:))]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (0:nanmax(numExchMets_3m(:)))+0.5; cbar.TickLabels = 0:nanmax(numExchMets_3m(:));
cbar.Label.String = 'Number of Exchanged Metabolites'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S3h', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S3h','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S4: What are the 9 transporters kept?

% Load Data
vars = {'min_trsptRxns*','biomass*','intlCon*'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_1m =  [  0,114,189]./256;
cmap_2m1 = [217, 83, 25]./256;
cmap_2m2 = [237,178, 32]./256;

% S4: Min Nt Reactions
[~,idx] = find([min_trsptRxns1m; min_trsptRxns2m1; min_trsptRxns2m2]); idx = unique(idx);
N1 = numel(find(biomass1(1,:)));
N2 = numel(find(biomass2(1,:)));
n = 17;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S4'; ax = axes(fig); clear h*
hold(ax,'on');
for metNum = 2:2:numel(idx)
    patch(ax,[0.5,1.5,1.5,0.5]+metNum-1,[0,0,100,100],[1,1,1].*0.75, 'FaceAlpha',0.5, 'EdgeColor','none');
end
h = bar(ax,100.*[min_trsptRxns1m(idx)./N1; min_trsptRxns2m1(idx)./N2; min_trsptRxns2m2(idx)./N2]',1); hold(ax,'off')
h(1).FaceColor = cmap_1m; h(2).FaceColor = cmap_2m1; h(3).FaceColor = cmap_2m2;
box(ax,'on'); grid(ax,'on'); axis(ax,'square'); xlim(ax,[0,numel(idx)]+0.5);
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
ax.XTick = 1:numel(idx); ax.XTickLabel = min_trsptRxns(idx); ax.XTickLabelRotation = 90;
legend(h,{'1 Strain','2 Strains (A)','2 Strains (B)'}, 'Location','NorthEastOutside', 'FontSize',legendSize, 'Box','Off'); clear h*
ax.Position = axPos;
xlabel(ax,'Transport Reactions', 'FontSize',xyLabelSize)
ylabel(ax,'Percentage of Simulations Kept When T_{TR}=9', 'FontSize',xyLabelSize)
title(ax,'Figure S4', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S4','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S5: Distance between 2-Strain Subnetworks

% Load Data
vars = {'jaccDist_2m','eucDist_2m','intlCon2','trsptCon2'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_Nt = flip(cbrewer('div','Spectral',numel(trsptCon2)-1));

% S5a: Jaccard Distance
n = 18;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S5a'; ax = axes(fig); clear h*
hold(ax,'on');
for ii = 1:numel(trsptCon2)-1
    h(ii) = plot(ax,intlCon2(1:end-2),jaccDist_2m(ii,1:end-2),'-', 'LineWidth',lineWidth, 'Color',cmap_Nt(ii,:));
end
hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,'square'); xlim(ax,intlLim); ylim(ax,[0, nanmax(jaccDist_2m(:))]);
cbar = colorbar(ax); caxis(ax,[0,numel(trsptCon2)-1]); colormap(ax,cmap_Nt);
cbar.Ticks = (1:(numel(trsptCon2)-1))-0.5; cbar.TickLabels = trsptCon2(1:end-1); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'T_{TR}'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
xlabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
ylabel(ax,'Jaccard Distance', 'FontSize',xyLabelSize)
title(ax,'Figure S5a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S5a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S5b: Euclidean Distance
n = 19;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S5b'; ax = axes(fig); clear h*
hold(ax,'on');
for ii = 1:numel(trsptCon2)-1
    h(ii) = plot(ax,intlCon2(1:end-2),eucDist_2m(ii,1:end-2),'-', 'LineWidth',lineWidth, 'Color',cmap_Nt(ii,:));
end
hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,'square'); xlim(ax,intlLim); ylim(ax,[0, nanmax(eucDist_2m(:))]);
cbar = colorbar(ax); caxis(ax,[0,numel(trsptCon2)-1]); colormap(ax,cmap_Nt);
cbar.Ticks = (1:(numel(trsptCon2)-1))-0.5; cbar.TickLabels = trsptCon2(1:end-1); cbar.FontSize = axesLabelSize;
cbar.Label.String = 'T_{TR}'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
xlabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
ylabel(ax,'Euclidean Distance', 'FontSize',xyLabelSize)
title(ax,'Figure S5b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S5b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S6: Number of Metabolites vs Jaccard Distance

% Load Data
vars = {'jaccDist*','numExchMets*'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_barPlot = 0.45.*ones(1,3);

% Calculate R^2 and Spearman Coeff
% 2-Strain Subnetworks
X = jaccDist_2m(:); Y = numExchMets_2m(:);
rm_idx = intersect(find(isnan(X)),find(isnan(Y))); % remove NaNs
X(rm_idx) = []; Y(rm_idx) = [];
rho2 = corr(X,Y, 'type','Pearson');
spear2 = corr(X,Y, 'type','Spearman'); clear X Y
% 3-Strain Subnetworks
X = jaccDist_3m(:); Y = numExchMets_3m(:);
rm_idx = intersect(find(isnan(X)),find(isnan(Y))); % remove NaNs
X(rm_idx) = []; Y(rm_idx) = [];
rho3 = corr(X,Y, 'type','Pearson');
spear3 = corr(X,Y, 'type','Spearman'); clear X Y

% S6a: Number of Exchanged Metabolites vs Jaccard Distance for 2-Strain Subnetworks
n = 20;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S6a'; ax = axes(fig); clear h*
scatter(ax,jaccDist_2m(:),numExchMets_2m(:),'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_barPlot, 'SizeData',50)
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,'tight','square'); xlim(ax,[0, nanmax([jaccDist_2m(:); jaccDist_3m(:)])]); ylim(ax,[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])]);
ax.Position = axPos; cbar.Position = cbarPos;
xlabel(ax,'Jaccard Distance', 'FontSize',xyLabelSize)
ylabel(ax,'Number of Exchanged Metabolites', 'FontSize',xyLabelSize)
title(ax,'Figure S6a', 'FontSize',titleSize)
% R^2
ah = annotation(fig,'textbox', 'String',['R^2=' num2str(round(rho2,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])-1, ah.Position(3:4)])
% Spearman Correlation
ah = annotation(fig,'textbox', 'String',['Spearman Coeff=' num2str(round(spear2,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])-3, ah.Position(3:4)])
print('DOLMN_Plots_Supp/Figure_S6a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S6b: Number of Exchanged Metabolites vs Jaccard Distance for 3-Strain Subnetworks
n = 21;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S6b'; ax = axes(fig); clear h*
scatter(ax,jaccDist_3m(:),numExchMets_3m(:),'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_barPlot, 'SizeData',50)
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,'tight','square'); xlim(ax,[0, nanmax([jaccDist_2m(:); jaccDist_3m(:)])]); ylim(ax,[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])]);
ax.Position = axPos; cbar.Position = cbarPos;
xlabel(ax,'Jaccard Distance', 'FontSize',xyLabelSize)
ylabel(ax,'Number of Exchanged Metabolites', 'FontSize',xyLabelSize)
title(ax,'Figure S6b', 'FontSize',titleSize)
% R^2
ah = annotation(fig,'textbox', 'String',['R^2=' num2str(round(rho3,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])-1, ah.Position(3:4)])
% Spearman Correlation
ah = annotation(fig,'textbox', 'String',['Spearman Coeff=' num2str(round(spear3,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[0, nanmax([numExchMets_2m(:); numExchMets_3m(:)])-3, ah.Position(3:4)])
print('DOLMN_Plots_Supp/Figure_S6b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars* rho* spear*

%% Figure S7: Predicting Metabolite Exchange

% Load Data
vars = {'exchMetsNum*','metNames','biomass*'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_barPlot = 0.45.*ones(1,3);
cmap_2m = [230, 85, 13]./256;
cmap_3m = [117,107,177]./256;

% Percentage of Simulations Exchanged
N1 = 100.*exchMetsNum_1m./numel(find(biomass1 ~= 0));
N2 = 100.*exchMetsNum_2m./numel(find(biomass2 ~= 0));
N3 = 100.*exchMetsNum_3m./numel(find(biomass3 ~= 0));

% Calculations
load('Ecoli_iJR904_newS.mat')
sol_wt = FBA(Ecoli,[],1);
intlMetsIdx = 1:numel(Ecoli.metNames); intlMetsIdx(identifyExchMets(Ecoli)) = [];
[~,bio_idx,~] = intersect(Ecoli.rxns,'BIOMASS_Ecoli');
cost = zeros(size(metNames));
numRxns = zeros(size(metNames));
bioCoeff = zeros(size(metNames));
for metNum = 1:numel(metNames)
    [~,rxn_idx,~] = intersect(erase(Ecoli.rxnNames,' exchange'), metNames{metNum});
    % Cost to Produce 1 mmol of Metabolite
    model = Ecoli;
    model.lb(rxn_idx) = 1;
    sol = FBA(model,[],1);
    cost(metNum) = 1 - sol.objectiveValue./sol_wt.objectiveValue;
    % Number of Reactions Associated with Metabolite
    [~,met_idx,~] = intersect(Ecoli.metNames(intlMetsIdx),metNames{metNum});
    met_idx = intlMetsIdx(met_idx);  
    numRxns(metNum) = numel(find(Ecoli.S(met_idx,:)));
    % Biomass Coefficient
    bioCoeff(metNum) = Ecoli.S(met_idx,bio_idx);
    clear model sol rxn_idx met_idx    
end
clear Ecoli sol* intlMetsIdx

% Spearman Coeff
spear2 = corr(N1,N2, 'type','Spearman');
spear3 = corr(N1,N3, 'type','Spearman');
% S7a: Percentage Exchanged vs Percentage Secreted Mets
n = 22;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S7a'; ax = axes(fig); clear h*
h(1) = scatter(ax,N1,N2,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_2m, 'SizeData',50); hold(ax,'on');
h(2) = scatter(ax,N1,N3,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_3m, 'SizeData',50); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0, 100, 0, 100],'square');
ax.Position = axPos;
legend(h,{['2 Strains, Spearman Coeff=' num2str(round(spear2,2))],['3 Strains, Spearman Coeff=' num2str(round(spear3,2))]}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'Secreted in 1-Strain Simulations (%)', 'FontSize',xyLabelSize)
ylabel(ax,'Exchanged in 2- & 3-Strain Simulations (%)', 'FontSize',xyLabelSize)
title(ax,'Figure S7a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S7a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])
clear spear*

% Spearman Coeff
rho23 = corr(N2,N3, 'type','Pearson');
spear23 = corr(N2,N3, 'type','Spearman');
% S7b: Percentage Exchanged Mets in 3 vs 2
n = 23;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S7b'; ax = axes(fig); clear h*
scatter(ax,N2,N3,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_barPlot, 'SizeData',50)
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0, 100, 0, 100],'square');
ax.Position = axPos;
xlabel(ax,'Exchanged in 2-Strain Simulations (%)', 'FontSize',xyLabelSize)
ylabel(ax,'Exchanged in 3-Strain Simulations (%)', 'FontSize',xyLabelSize)
title(ax,'Figure S7b', 'FontSize',titleSize)
% R^2
ah = annotation(fig,'textbox', 'String',['R^2=' num2str(round(rho23,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[70, 6, ah.Position(3:4)])
% Spearman Correlation
ah = annotation(fig,'textbox', 'String',['Spearman Coeff=' num2str(round(spear23,2))], 'FontSize',legendSize, 'LineStyle','none');
set(ah, 'Parent',ax, 'Position',[70, 10, ah.Position(3:4)])
print('DOLMN_Plots_Supp/Figure_S7b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])
clear rho* spear*

% Spearman Coeff
spear2 = corr(cost,N2, 'type','Spearman');
spear3 = corr(cost,N3, 'type','Spearman');
% S7c: Number of Exchanged Mets vs Cost to Produce
n = 24;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S7c'; ax = axes(fig); clear h*
h(1) = scatter(ax,cost,N2,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_2m, 'SizeData',50); hold(ax,'on');
h(2) = scatter(ax,cost,N3,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_3m, 'SizeData',50); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0, max(cost), 0, 100],'square');
ax.Position = axPos;
legend(h,{['2 Strains, Spearman Coeff=' num2str(round(spear2,2))],['3 Strains, Spearman Coeff=' num2str(round(spear3,2))]}, 'Location','NorthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'Cost to Produce 1 mmol of Metabolite', 'FontSize',xyLabelSize)
ylabel(ax,'Exchanged in 2- & 3-Strain Simulations (%)', 'FontSize',xyLabelSize)
title(ax,'Figure S7c', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S7c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])
clear spear*

% Spearman Coeff
spear2 = corr(numRxns,N2, 'type','Spearman');
spear3 = corr(numRxns,N3, 'type','Spearman');
% S7d: Number of Exchanged Mets vs Number of Reactions
n = 25;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S7d'; ax = axes(fig); clear h*
h(1) = scatter(ax,numRxns,N2,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_2m, 'SizeData',50); hold(ax,'on');
h(2) = scatter(ax,numRxns,N3,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_3m, 'SizeData',50); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[0, max(numRxns), 0, 100],'square');
ax.Position = axPos;
legend(h,{['2 Strains, Spearman Coeff=' num2str(round(spear2,2))],['3 Strains, Spearman Coeff=' num2str(round(spear3,2))]}, 'Location','NorthEast', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'Number of Reactions Associated with Metabolite', 'FontSize',xyLabelSize)
ylabel(ax,'Exchanged in 2- & 3-Strain Simulations (%)', 'FontSize',xyLabelSize)
title(ax,'Figure S7d', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S7d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])
clear spear*

% R2
rho2 = corr(bioCoeff,N2, 'type','Pearson');
rho3 = corr(bioCoeff,N3, 'type','Pearson');
% Spearman Coeff
spear2 = corr(bioCoeff,N2, 'type','Spearman');
spear3 = corr(bioCoeff,N3, 'type','Spearman');
% S7e: Number of Exchanged Mets vs Biomass Coefficient
n = 26;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S7e'; ax = axes(fig); clear h*
h(1) = scatter(ax,bioCoeff,N2,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_2m, 'SizeData',50); hold(ax,'on');
h(2) = scatter(ax,bioCoeff,N3,'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmap_3m, 'SizeData',50); hold(ax,'off');
box(ax,'on'); grid(ax,'on'); ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
axis(ax,[min(bioCoeff), 0, 0, 100],'square');
ax.Position = axPos;
legend(h,{['2 Strains, Spearman Coeff=' num2str(round(spear2,2)) ', R^2=' num2str(round(rho2,2))],['3 Strains, Spearman Coeff=' num2str(round(spear3,2)) ', R^2=' num2str(round(rho3,2))]}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'Biomass Coefficient', 'FontSize',xyLabelSize)
ylabel(ax,'Exchanged in 2- & 3-Strain Simulations (%)', 'FontSize',xyLabelSize)
title(ax,'Figure S7e', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S7e','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])
clear rho* spear*

% Clear Variables
clear(vars{:})
clear cmap* vars* N*
clear cost numRxns bioCoeff

%% Figure S8: PCA

% Load Data
vars = {'trsptCon2','score2m1','score2m2','explained','Nt2m','Ni2m'}; % PCA
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% S8: PCA
n = 27;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S8'; clear h*
for trsptNum = 1:numel(trsptCon2)-1
    idx_t = find(Nt2m == trsptCon2(trsptNum));
    idx_t(Ni2m(idx_t) > 285) = []; % remove Ni > 285
    cmap_pca = flip(cbrewer('div','Spectral',numel(idx_t)));
    
    ax(trsptNum) = subplot(5,8,trsptNum);
    plot(ax(trsptNum),score2m1(idx_t,1),score2m1(idx_t,2),'k-', 'LineWidth',0.5*lineWidth); hold(ax(trsptNum),'on');
    plot(ax(trsptNum),score2m2(idx_t,1),score2m2(idx_t,2),'k--', 'LineWidth',0.5*lineWidth);
    for ii = 1:numel(idx_t) % increasing Ni
        plot(ax(trsptNum),score2m1(idx_t(ii),1),score2m1(idx_t(ii),2),'v', 'MarkerSize',5, 'MarkerFaceColor',cmap_pca(ii,:), 'MarkerEdgeColor','none')
        plot(ax(trsptNum),score2m2(idx_t(ii),1),score2m2(idx_t(ii),2),'^', 'MarkerSize',5, 'MarkerFaceColor',cmap_pca(ii,:), 'MarkerEdgeColor','none')
    end
    hold(ax(trsptNum),'off');
    box(ax(trsptNum),'on'); grid(ax(trsptNum),'on'); ax(trsptNum).FontSize = 0.5*axesLabelSize; ax(trsptNum).TickLength = tickLength.*[1 2];
    title(ax(trsptNum),int2str(trsptCon2(trsptNum)), 'FontSize',0.5*titleSize)
    
    clear idx_t
end
% Colorbar
temp = subplot(5,8,trsptNum+1); temp.Visible = 'off';
cbar = colorbar(temp); colormap(temp,cmap_pca);
cbar.XTick = cbar.XTick([1,end]); cbar.XTickLabel = {'Small T_{IN}'; 'Large T_{IN}'};
linkaxes(ax,'xy'); xlim(ax(trsptNum),[-200,300]); ylim(ax(trsptNum),[-50,250])
[~,hx]=suplabel(['PC1: ' num2str(round(explained(1),1)) '%'],'x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel(['PC2: ' num2str(round(explained(2),1)) '%'],'y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S8','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S8','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S9-10: Metabolic Differentiation

% Load Data
vars = {'trsptCon','intlCon','jaccDist_2m','pathwayEucDist','pathwayJaccDist','pathwayNames','metNames','exchMets_2m','numExchMets_2m','x1','y1','x2','y2'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_interp.mat']),vars{:})

% Colormaps
cmap_jaccDist = cbrewer('seq','Blues',256);
cmap_eucDist = cbrewer('seq','Reds',256);

climits = [0, max(arrayfun(@(x) nanmax(pathwayJaccDist{x}(:)),1:numel(pathwayNames)))];
% S9: Jaccard Distance
n = 28;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S9'; clear h*
for pathNum = 1:numel(pathwayNames)
    ax = subplot(5,6,pathNum);
    patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
    g = imagesc(ax,trsptCon,intlCon,pathwayJaccDist{pathNum}'); g.AlphaData = ~isnan(pathwayJaccDist{pathNum}'); clear g
    plot(ax,x1,y1,'k--', 'LineWidth',0.33*lineWidth);
    plot(ax,x2,y2,'k-', 'LineWidth',0.33*lineWidth); hold(ax,'off');
    box(ax,'on'); grid(ax,'on'); ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    axis(ax,[trsptLim, intlLim],'square');
    colormap(ax,cmap_jaccDist); caxis(ax,climits);
    title(ax,pathwayNames{pathNum}, 'FontSize',0.5*titleSize)
end
% Colorbar
temp = subplot(5,6,pathNum+1); temp.Visible = 'off';
cbar = colorbar(temp); colormap(temp,cmap_jaccDist); cbar.FontSize = 0.5*axesLabelSize;
cbar.Label.String = 'Jaccard Distance'; cbar.Label.FontSize = 0.5*xyLabelSize;
[~,hx]=suplabel('T_{TR}','x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel('T_{IN}','y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S9','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S9','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S10: Euclidean Distance
n = 29;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S10'; clear h*
for pathNum = 1:numel(pathwayNames)
    ax = subplot(5,6,pathNum);
    patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
    g = imagesc(ax,trsptCon,intlCon,pathwayEucDist{pathNum}'); g.AlphaData = ~isnan(pathwayEucDist{pathNum}'); clear g
    plot(ax,x1,y1,'k--', 'LineWidth',0.33*lineWidth);
    plot(ax,x2,y2,'k-', 'LineWidth',0.33*lineWidth); hold(ax,'off');
    box(ax,'on'); grid(ax,'on'); ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    axis(ax,[trsptLim, intlLim],'square');
    colormap(ax,cmap_eucDist);
    if nanmax(pathwayEucDist{pathNum}(:)) > 1E-3
        caxis(ax,[0,nanmax(pathwayEucDist{pathNum}(:))]);
        cbar = colorbar(ax); cbar.FontSize = 0.5*axesLabelSize;
        if nanmax(pathwayEucDist{pathNum}(:)) > 1
            cbar.Ticks = caxis; cbar.XTickLabel = cellstr(int2str(round(caxis')));
        else
            cbar.Ticks = caxis; cbar.XTickLabel = cellstr(num2str(round(caxis',2)));
        end
    else
        caxis(ax,[0,100]);
        cbar = colorbar(ax); cbar.FontSize = 0.5*axesLabelSize;
        cbar.Ticks = caxis; cbar.XTickLabel = {'0';'0'};
    end
    title(ax,pathwayNames{pathNum}, 'FontSize',0.5*titleSize)
end
[~,hx]=suplabel('T_{TR}','x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel('T_{IN}','y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S10','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S10','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*

%% Figure S11-13: Metabolite secretion/exchange profiles

% Load Data
vars1 = {'trsptCon','intlCon','metNames','exchMets*','x1','y1','x2','y2','x3','y3'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_interp.mat']),vars1{:})
vars2 = {'P1_phi','P2_phi','P3_phi','phi*','exchMetsNum*','exchMetNames'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars2{:})

% Colormaps
cmap_mets = 0.45.*ones(1,3);

% Re-Order by Hierarchical Clustering - 1 Strain
phi = phi1(P1_phi,P1_phi);
tempNames = exchMetNames(P1_phi);
rm_idx = find(ismember(phi,zeros(1,size(phi,1)),'rows')); % remove all zeros
tempNames(rm_idx) = [];
[~,~,idx_met1m] = intersect(tempNames,metNames, 'stable'); idx_met1m = flip(idx_met1m);
clear phi tempNames

% Re-Order by Hierarchical Clustering - 2 Strains
phi = phi2(P2_phi,P2_phi);
tempNames = exchMetNames(P2_phi);
rm_idx = find(ismember(phi,zeros(1,size(phi,1)),'rows')); % remove all zeros
tempNames(rm_idx) = [];
[~,~,idx_met2m] = intersect(tempNames,metNames, 'stable'); idx_met2m = flip(idx_met2m);
clear phi tempNames

% Re-Order by Hierarchical Clustering - 3 Strains
phi = phi3(P3_phi,P3_phi);
tempNames = exchMetNames(P3_phi);
rm_idx = find(ismember(phi,zeros(1,size(phi,1)),'rows')); % remove all zeros
tempNames(rm_idx) = [];
[~,~,idx_met3m] = intersect(tempNames,metNames, 'stable'); idx_met3m = flip(idx_met3m);
clear phi tempNames

% S11: 1-Strain Subnetworks
n = 30;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S11'; clear h*
for metNum = 1:numel(idx_met1m)
    ax = subplot(3,5,metNum);
    patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
    g = imagesc(ax,trsptCon,intlCon,exchMets_1m{idx_met1m(metNum)}'); g.AlphaData = ~isnan(exchMets_1m{idx_met1m(metNum)}'); clear g
    plot(ax,x1,y1,'k--', 'LineWidth',0.33*lineWidth); hold(ax,'off');
    box(ax,'on'); grid(ax,'on'); ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    axis(ax,[trsptLim, intlLim],'square');
    colormap(ax,cmap_mets); caxis(ax,[0, 1]);
    title(ax,metNames{idx_met1m(metNum)}, 'FontSize',0.5*titleSize)
end
[~,hx]=suplabel('T_{TR}','x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel('T_{IN}','y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S11','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S11','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S12: 2-Strain Subnetworks
n = 31;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S12'; clear h*
for metNum = 1:numel(idx_met2m)
    ax = subplot(5,8,metNum);
    patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
    g = imagesc(ax,trsptCon,intlCon,exchMets_2m{idx_met2m(metNum)}'); g.AlphaData = ~isnan(exchMets_2m{idx_met2m(metNum)}'); clear g
    plot(ax,x1,y1,'k--', 'LineWidth',0.33*lineWidth);
    plot(ax,x2,y2,'k-', 'LineWidth',0.33*lineWidth); hold(ax,'off');
    box(ax,'on'); grid(ax,'on'); ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    axis(ax,[trsptLim, intlLim],'square');
    colormap(ax,cmap_mets); caxis(ax,[0, 1]);
    title(ax,metNames{idx_met2m(metNum)}, 'FontSize',0.5*titleSize)
end
[~,hx]=suplabel('T_{TR}','x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel('T_{IN}','y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S12','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S12','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S13: 3-Strain Subnetworks
n = 32;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S13'; clear h*
for metNum = 1:numel(idx_met3m)
    ax = subplot(5,8,metNum);
    patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on'); 
    g = imagesc(ax,trsptCon,intlCon,exchMets_3m{idx_met3m(metNum)}'); g.AlphaData = ~isnan(exchMets_3m{idx_met3m(metNum)}'); clear g
    plot(ax,x1,y1,'k--', 'LineWidth',0.33*lineWidth);
    plot(ax,x2,y2,'k-', 'LineWidth',0.33*lineWidth);
    plot(ax,x3,y3,'k:', 'LineWidth',0.33*lineWidth); hold(ax,'off');
    box(ax,'on'); grid(ax,'on'); ax.FontSize = 0.5*axesLabelSize; ax.TickLength = tickLength.*[1 2];
    axis(ax,[trsptLim, intlLim],'square');
    colormap(ax,cmap_mets); caxis(ax,[0, 1]);
    title(ax,metNames{idx_met3m(metNum)}, 'FontSize',0.5*titleSize)
end
[~,hx]=suplabel('T_{TR}','x'); hx.FontSize = xyLabelSize; hx.Position = [0.5,0,0];
[~,hy]=suplabel('T_{IN}','y'); hy.FontSize = xyLabelSize; hy.Position = [0,0.5,0];
[~,ht]=suplabel('Figure S13','t'); ht.FontSize = titleSize; ht.Position = [0.5,1.006,0.5];
print('DOLMN_Plots_Supp/Figure_S13','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars1{:})
clear(vars2{:})
clear cmap* vars*
clear metNames* idx* rho

%% Figure S14: C, N, O

% Load Data
vars = {'trsptCon','intlCon','C_*','N_*','O_*','x1','y1','x2','y2','x3','y3'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_interp.mat']),vars{:})

% Colormaps
cmap_1m = [[117,107,177]; [240,240,240]]./256;
cmap_2m = [[230,85,13];              [240,240,240]]./256;
cmap_3m = [[230,85,13]; [49,163,84]; [240,240,240]]./256;

% Colormap Labels
clabel_1m = {'Does Not Use'; 'Uses'};
clabel_2m = {'1 Strain Uses';                  'All Strains Use'};
clabel_3m = {'1 Strain Uses'; '2 Strains Use'; 'All Strains Use'};

% S14a: C2
n = 33;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14a'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,C_2m'); g.AlphaData = ~isnan(C_2m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_2m); caxis(ax,[1, 3]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:2)+0.5; cbar.TickLabels = clabel_2m;
cbar.Label.String = 'Glucose'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S14b: C3
n = 34;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14b'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,C_3m'); g.AlphaData = ~isnan(C_3m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_3m); caxis(ax,[1, 4]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:3)+0.5; cbar.TickLabels = clabel_3m;
cbar.Label.String = 'Glucose'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S14c: N2
n = 35;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14a'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,N_2m'); g.AlphaData = ~isnan(N_2m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_2m); caxis(ax,[1, 3]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:2)+0.5; cbar.TickLabels = clabel_2m;
cbar.Label.String = 'Ammonia'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14c', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14c','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S14d: N3
n = 36;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14d'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,N_3m'); g.AlphaData = ~isnan(N_3m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_3m); caxis(ax,[1, 4]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:3)+0.5; cbar.TickLabels = clabel_3m;
cbar.Label.String = 'Ammonia'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14d', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14d','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S14e: O2
n = 37;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14e'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,O_2m'); g.AlphaData = ~isnan(O_2m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_2m); caxis(ax,[1, 3]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:2)+0.5; cbar.TickLabels = clabel_2m;
cbar.Label.String = 'Oxygen'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains'}, 'Location','SouthEast', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14e', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14e','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S14f: O3
n = 38;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S14f'; ax = axes(fig); clear h*
patch(ax,[0,8.5,8.5,0],[0,0,285,285],[1,1,1].*0.75, 'FaceAlpha',0.75, 'EdgeColor','none'); hold(ax,'on');
g = imagesc(ax,trsptCon,intlCon,O_3m'); g.AlphaData = ~isnan(O_3m'); clear g;
h(1) = plot(ax,x1,y1,'k--', 'LineWidth',lineWidth);
h(2) = plot(ax,x2,y2,'k-', 'LineWidth',lineWidth);
h(3) = plot(ax,x3,y3,'k:', 'LineWidth',lineWidth); hold(ax,'off');
axis(ax,[trsptLim, intlLim],'square'); box(ax,'on'); grid(ax,'on');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2]; ax.YDir = 'normal';
colormap(ax,cmap_3m); caxis(ax,[1, 4]);
cbar = colorbar(ax); cbar.FontSize = axesLabelSize;
cbar.Ticks = (1:3)+0.5; cbar.TickLabels = clabel_3m;
cbar.Label.String = 'Oxygen'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
legend(h,{'1 Strain','2 Strains','3 Strains'}, 'Location','SouthWest', 'FontSize',legendSize, 'Box','Off'); clear h*
xlabel(ax,'T_{TR}', 'FontSize',xyLabelSize)
ylabel(ax,'T_{IN}', 'FontSize',xyLabelSize)
title(ax,'Figure S14f', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S14f','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*
clear clabel*

%% Figure S15: Phi Coefficient

% Load Data
vars = {'phi2','P2_rho','exchMetsNum_2m','exchMetNames','mediumMets','metNames','biomass2'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_phi = flip(cbrewer('div','RdBu',256));

% Hierarchical Clustering
phi = phi2(P2_rho,P2_rho);
exchMetNames = exchMetNames(P2_rho);
rm_idx = find(ismember(phi,zeros(1,size(phi,1)),'rows')); % remove all zeros
phi(rm_idx,:) = []; phi(:,rm_idx) = []; exchMetNames(rm_idx) = [];
[~,rm_idx,~] = intersect(exchMetNames,mediumMets); % remove medium metabolites
phi(rm_idx,:) = []; phi(:,rm_idx) = []; exchMetNames(rm_idx) = [];

% S15: Phi Coefficient
n = 39;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S15'; ax = axes(fig); clear h*
imagesc(ax,phi);
box(ax,'on'); grid(ax,'off'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
cbar = colorbar(ax); caxis(ax,[-1,1]); colormap(ax,cmap_phi);
cbar.Label.String = 'Phi Coefficient'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
ax.XTickLabelRotation = 90; ax.TickLabelInterpreter = 'none';
ax.XTick = 1:numel(exchMetNames); ax.XTickLabel = exchMetNames;
ax.YTick = 1:numel(exchMetNames); ax.YTickLabel = exchMetNames;
title(ax,'Figure S15', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S15','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*
clear phi

%% Figure S18: Point-Biserial Correlation Between Jaccard Distance and Exchanged Metabolite Landscape

% Load Data
vars = {'corr_pathMet*','P_path*','pathwayNames','P_mets*','metNames'};
load(fullfile('DOLMN_Parsed',[loadDataName_full '_summary.mat']),vars{:})

% Colormaps
cmap_corr = flip(cbrewer('div','RdBu',256));

% Hierarchical Clustering
corrA = corr_pathMet_euc(P_pathEuc,P_metsEuc);
corrB = corr_pathMet_jacc(P_pathEuc,P_metsEuc);

% S18a: Point-Biserial Correlation (Euclidean)
n = 40;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S18a'; ax = axes(fig); clear h*
imagesc(ax,corrA);
box(ax,'on'); grid(ax,'off'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
cbar = colorbar(ax); caxis(ax,[-1,1]); colormap(ax,cmap_corr);
cbar.Label.String = 'Point-Biserial Correlation'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
ax.XTickLabelRotation = 90; ax.TickLabelInterpreter = 'none';
ax.XTick = 1:numel(metNames); ax.XTickLabel = metNames(P_metsEuc);
ax.YTick = 1:numel(pathwayNames); ax.YTickLabel = pathwayNames(P_pathEuc);
xlabel(ax,'Metabolites', 'FontSize',xyLabelSize)
ylabel(ax,'Pathways', 'FontSize',xyLabelSize)
title(ax,'Figure S18a', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S18a','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% S18a: Point-Biserial Correlation (Jaccard)
n = 41;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'S18b'; ax = axes(fig); clear h*
imagesc(ax,corrB);
box(ax,'on'); grid(ax,'off'); axis(ax,'square');
ax.FontSize = axesLabelSize; ax.TickLength = tickLength.*[1 2];
cbar = colorbar(ax); caxis(ax,[-1,1]); colormap(ax,cmap_corr);
cbar.Label.String = 'Point-Biserial Correlation'; cbar.Label.FontSize = xyLabelSize;
ax.Position = axPos; cbar.Position = cbarPos;
ax.XTickLabelRotation = 90; ax.TickLabelInterpreter = 'none';
ax.XTick = 1:numel(metNames); ax.XTickLabel = metNames(P_metsEuc);
ax.YTick = 1:numel(pathwayNames); ax.YTickLabel = pathwayNames(P_pathEuc);
xlabel(ax,'Metabolites', 'FontSize',xyLabelSize)
ylabel(ax,'Pathways', 'FontSize',xyLabelSize)
title(ax,'Figure S18b', 'FontSize',titleSize)
print('DOLMN_Plots_Supp/Figure_S18b','-dpdf','-bestfit','-r300','-painters',['-f' int2str(fig.Number)])

% Clear Variables
clear(vars{:})
clear cmap* vars*
clear corr*










