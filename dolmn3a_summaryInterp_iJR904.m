clear variables
close all
clc

%% Load Data

loadDataName = 'Ecoli_iJR904';

% Load Metabolic Model
load('Ecoli_iJR904_newS.mat')

% Load DOLMN Data
load(['DOLMN_Parsed\iJR904\' loadDataName '.mat'])

%% Interpolate Biomass

trsptCon = min([trsptCon1, trsptCon2, trsptCon3]):45;
intlCon = min([intlCon1, intlCon2, intlCon3]):285;
[intlCon_mat,trsptCon_mat] = meshgrid(intlCon,trsptCon);

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.biomass,1:numel(model1), 'Uni',false)');
biomass1 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.biomass,1:numel(model2), 'Uni',false)');
biomass2 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.biomass,1:numel(model3), 'Uni',false)');
biomass3 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);

%% Biomass Flux Growth Boundary

% Matrices
biomass1_nan = biomass1; biomass1_nan(biomass1_nan == 0) = NaN;
biomass2_nan = biomass2; biomass2_nan(biomass2_nan == 0) = NaN;
biomass3_nan = biomass3; biomass3_nan(biomass3_nan == 0) = NaN;

% 1 Model
B = bwboundaries(~isnan(biomass1_nan));
x1 = trsptCon(B{1}(:,1)); y1 = intlCon(B{1}(:,2));
rm_idx = find(x1 == max(x1) & y1 > min(y1)); x1(rm_idx) = []; y1(rm_idx) = [];
rm_idx = find(y1 == max(y1) & x1 > min(x1)); x1(rm_idx) = []; y1(rm_idx) = [];
idx = find(x1 == max(x1)); x1 = [x1(idx:end), x1(1:idx-1)]; y1 = [y1(idx:end), y1(1:idx-1)];
z1 = arrayfun(@(k) biomass1(trsptCon == x1(k),intlCon == y1(k)),1:numel(x1));
% 2 Models
B = bwboundaries(~isnan(biomass2_nan));
x2 = trsptCon(B{1}(:,1)); y2 = intlCon(B{1}(:,2));
rm_idx = find(x2 == max(x2) & y2 > min(y2)); x2(rm_idx) = []; y2(rm_idx) = [];
rm_idx = find(y2 == max(y2) & x2 > min(x2)); x2(rm_idx) = []; y2(rm_idx) = [];
idx = find(x2 == max(x2)); x2 = [x2(idx:end), x2(1:idx-1)]; y2 = [y2(idx:end), y2(1:idx-1)];
z2 = arrayfun(@(k) biomass2(trsptCon == x2(k),intlCon == y2(k)),1:numel(x2));
% 3 Models
B = bwboundaries(~isnan(biomass3_nan));
x3 = trsptCon(B{1}(:,1)); y3 = intlCon(B{1}(:,2));
rm_idx = find(x3 == max(x3) & y3 > min(y3)); x3(rm_idx) = []; y3(rm_idx) = [];
rm_idx = find(y3 == max(y3) & x3 > min(x3)); x3(rm_idx) = []; y3(rm_idx) = [];
idx = find(x3 == max(x3)); x3 = [x3(idx:end), x3(1:idx-1)]; y3 = [y3(idx:end), y3(1:idx-1)];
z3 = arrayfun(@(k) biomass3(trsptCon == x3(k),intlCon == y3(k)),1:numel(x3));

x1(2:end) = x1(2:end) - 0.5; x2(2:end) = x2(2:end) - 0.5; x3(2:end) = x3(2:end) - 0.5;
y1(1:end-1) = y1(1:end-1) - 0.5; y2(1:end-1) = y2(1:end-1) - 0.5; y3(1:end-1) = y3(1:end-1) - 0.5;

%% Interpolate Exchanged/Secreted Metabolites

% 1 Model
exchMetsNum_1m = zeros(Ecoli.N_u,1);
exchMets_1m = cell(Ecoli.N_u,1);
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) sum(model1{x}{1}.flux(model1{x}{1}.exch_idx,:) > 0),1:numel(model1), 'Uni',false)');
numExchMets_1m = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); numExchMets_1m(biomass1 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(model1{x}{1}.exch_idx(metNum),:),1:numel(model1), 'Uni',false)');
    V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V < 0) = 0;
    exchMets_1m{metNum} = double(logical(V));
    exchMetsNum_1m(metNum) = sum(exchMets_1m{metNum}(:));
    exchMets_1m{metNum}(exchMets_1m{metNum} == 0) = NaN;
    exchMets_1m{metNum}(biomass1 == 0) = NaN;
end

% 2 Models
exchMetsNum_2m = zeros(Ecoli.N_u,1);
exchMets_2m = cell(Ecoli.N_u,1);
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) sum(model2{x}{1}.exchMets.model1_to_model2 + model2{x}{1}.exchMets.model2_to_model1),1:numel(model2), 'Uni',false)');
numExchMets_2m = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); numExchMets_2m(biomass2 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    V = cell2mat(arrayfun(@(x) model2{x}{1}.exchMets.model1_to_model2(model2{x}{1}.exch_idx(metNum),:) + model2{x}{1}.exchMets.model2_to_model1(model2{x}{1}.exch_idx(metNum),:),1:numel(model2), 'Uni',false)');
    V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
    exchMets_2m{metNum} = double(logical(V));
    exchMetsNum_2m(metNum) = sum(exchMets_2m{metNum}(:));
    exchMets_2m{metNum}(exchMets_2m{metNum} == 0) = NaN;
    exchMets_2m{metNum}(biomass2 == 0) = NaN;
end

% 3 Models
exchMetsNum_3m = zeros(Ecoli.N_u,1);
exchMets_3m = cell(Ecoli.N_u,1);
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) sum(logical(model3{x}{1}.exchMets.model1_to_model2 + model3{x}{1}.exchMets.model2_to_model1 + ...
    model3{x}{1}.exchMets.model1_to_model3 + model3{x}{1}.exchMets.model3_to_model1 + ...
    model3{x}{1}.exchMets.model2_to_model3 + model3{x}{1}.exchMets.model3_to_model2)),1:numel(model3), 'Uni',false)');
numExchMets_3m = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); numExchMets_3m(biomass3 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    V = cell2mat(arrayfun(@(x) model3{x}{1}.exchMets.model1_to_model2(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model2_to_model1(model3{x}{1}.exch_idx(metNum),:) + ...
        model3{x}{1}.exchMets.model1_to_model3(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model3_to_model1(model3{x}{1}.exch_idx(metNum),:) + ...
        model3{x}{1}.exchMets.model2_to_model3(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model3_to_model2(model3{x}{1}.exch_idx(metNum),:),1:numel(model3), 'Uni',false)');
    V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
    exchMets_3m{metNum} = double(logical(V));
    exchMetsNum_3m(metNum) = sum(exchMets_3m{metNum}(:));
    exchMets_3m{metNum}(exchMets_3m{metNum} == 0) = NaN;
    exchMets_3m{metNum}(biomass3 == 0) = NaN;
end

%% Interpolate %% Jaccard & Euclidean Distance

% 2 Models - Jaccard
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.int',model2{x}{2}.int','jaccard')),1:numel(model2), 'Uni',false))';
jaccDist_2m = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); jaccDist_2m(biomass2 == 0) = NaN;
% 2 Models - Euclidean
V = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.flux',model2{x}{2}.flux','euclidean')),1:numel(model2), 'Uni',false))';
eucDist_2m = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); eucDist_2m(biomass2 == 0) = NaN;

% 3 Models - Jaccard
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.int'),double(model3{x}{2}.int'),'jaccard')),1:numel(model3), 'Uni',false))';
jaccDist_3m12 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.int'),double(model3{x}{3}.int'),'jaccard')),1:numel(model3), 'Uni',false))';
jaccDist_3m13 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{2}.int'),double(model3{x}{3}.int'),'jaccard')),1:numel(model3), 'Uni',false))';
jaccDist_3m23 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
% average jaccard's distance
jaccDist_3m = (jaccDist_3m12 + jaccDist_3m13 + jaccDist_3m23)./3; jaccDist_3m(biomass3 == 0) = NaN;
jaccDist_3m12(biomass3 == 0) = NaN; jaccDist_3m13(biomass3 == 0) = NaN; jaccDist_3m23(biomass3 == 0) = NaN;
% 3 Models - Euclidean
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.flux'),double(model3{x}{2}.flux'),'euclidean')),1:numel(model3), 'Uni',false))';
eucDist_3m12 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.flux'),double(model3{x}{3}.flux'),'euclidean')),1:numel(model3), 'Uni',false))';
eucDist_3m13 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
V = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{2}.flux'),double(model3{x}{3}.flux'),'euclidean')),1:numel(model3), 'Uni',false))';
eucDist_3m23 = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0);
% average Euclidean distance
eucDist_3m = (eucDist_3m12 + eucDist_3m13 + eucDist_3m23)./3; eucDist_3m(biomass3 == 0) = NaN;
eucDist_3m12(biomass3 == 0) = NaN; eucDist_3m13(biomass3 == 0) = NaN; eucDist_3m23(biomass3 == 0) = NaN;

%% Interpolate Where Models Don't Use Glucose

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_glc__D_e');

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_1m = double(logical(V)); C_1m(biomass1 == 0) = NaN;

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_2m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_2m2 = double(logical(V));
C_2m = zeros(numel(trsptCon),numel(intlCon));
C_2m(C_2m1==1 & C_2m2==1) = 2;
C_2m(C_2m1==1 & C_2m2==0) = 1;
C_2m(C_2m1==0 & C_2m2==1) = 1;
C_2m(biomass2 == 0) = NaN;

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_3m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_3m2 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
C_3m3 = double(logical(V));
C_3m = zeros(numel(trsptCon),numel(intlCon));
C_3m(C_3m1==1 & C_3m2==1 & C_3m3==1) = 3; % all
C_3m(C_3m1==1 & C_3m2==0 & C_3m3==0) = 1; % 1 only
C_3m(C_3m1==0 & C_3m2==1 & C_3m3==0) = 1; % 2 only
C_3m(C_3m1==0 & C_3m2==0 & C_3m3==1) = 1; % 3 only
C_3m(C_3m1==1 & C_3m2==1 & C_3m3==0) = 2; % 1&2
C_3m(C_3m1==1 & C_3m2==0 & C_3m3==1) = 2; % 1&3
C_3m(C_3m1==0 & C_3m2==1 & C_3m3==1) = 2; % 2&3
C_3m(biomass3 == 0) = NaN;

%% Interpolate Where Models Don't Use Ammonia

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_nh4_e');

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_1m = double(logical(V)); N_1m(biomass1 == 0) = NaN;

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_2m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_2m2 = double(logical(V));
N_2m = zeros(numel(trsptCon),numel(intlCon));
N_2m(N_2m1==1 & N_2m2==1) = 2;
N_2m(N_2m1==1 & N_2m2==0) = 1;
N_2m(N_2m1==0 & N_2m2==1) = 1;
N_2m(biomass2 == 0) = NaN;

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_3m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_3m2 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
N_3m3 = double(logical(V));
N_3m = zeros(numel(trsptCon),numel(intlCon));
N_3m(N_3m1==1 & N_3m2==1 & N_3m3==1) = 3; % all
N_3m(N_3m1==1 & N_3m2==0 & N_3m3==0) = 1; % 1 only
N_3m(N_3m1==0 & N_3m2==1 & N_3m3==0) = 1; % 2 only
N_3m(N_3m1==0 & N_3m2==0 & N_3m3==1) = 1; % 3 only
N_3m(N_3m1==1 & N_3m2==1 & N_3m3==0) = 2; % 1&2
N_3m(N_3m1==1 & N_3m2==0 & N_3m3==1) = 2; % 1&3
N_3m(N_3m1==0 & N_3m2==1 & N_3m3==1) = 2; % 2&3
N_3m(biomass3 == 0) = NaN;

%% Interpolate Where Models Don't Use Phosphate

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_pi_e');

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_1m = double(logical(V)); P_1m(biomass1 == 0) = NaN;

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_2m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_2m2 = double(logical(V));
P_2m = zeros(numel(trsptCon),numel(intlCon));
P_2m(P_2m1==1 & P_2m2==1) = 2;
P_2m(P_2m1==1 & P_2m2==0) = 1;
P_2m(P_2m1==0 & P_2m2==1) = 1;
P_2m(biomass2 == 0) = NaN;

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_3m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_3m2 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
P_3m3 = double(logical(V));
P_3m = zeros(numel(trsptCon),numel(intlCon));
P_3m(P_3m1==1 & P_3m2==1 & P_3m3==1) = 3; % all
P_3m(P_3m1==1 & P_3m2==0 & P_3m3==0) = 1; % 1 only
P_3m(P_3m1==0 & P_3m2==1 & P_3m3==0) = 1; % 2 only
P_3m(P_3m1==0 & P_3m2==0 & P_3m3==1) = 1; % 3 only
P_3m(P_3m1==1 & P_3m2==1 & P_3m3==0) = 2; % 1&2
P_3m(P_3m1==1 & P_3m2==0 & P_3m3==1) = 2; % 1&3
P_3m(P_3m1==0 & P_3m2==1 & P_3m3==1) = 2; % 2&3
P_3m(biomass3 == 0) = NaN;

%% Interpolate Where Models Don't Use Sulfate

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_so4_e');

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_1m = double(logical(V)); S_1m(biomass1 == 0) = NaN;

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_2m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_2m2 = double(logical(V));
S_2m = zeros(numel(trsptCon),numel(intlCon));
S_2m(S_2m1==1 & S_2m2==1) = 2;
S_2m(S_2m1==1 & S_2m2==0) = 1;
S_2m(S_2m1==0 & S_2m2==1) = 1;
S_2m(biomass2 == 0) = NaN;

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_3m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_3m2 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
S_3m3 = double(logical(V));
S_3m = zeros(numel(trsptCon),numel(intlCon));
S_3m(S_3m1==1 & S_3m2==1 & S_3m3==1) = 3; % all
S_3m(S_3m1==1 & S_3m2==0 & S_3m3==0) = 1; % 1 only
S_3m(S_3m1==0 & S_3m2==1 & S_3m3==0) = 1; % 2 only
S_3m(S_3m1==0 & S_3m2==0 & S_3m3==1) = 1; % 3 only
S_3m(S_3m1==1 & S_3m2==1 & S_3m3==0) = 2; % 1&2
S_3m(S_3m1==1 & S_3m2==0 & S_3m3==1) = 2; % 1&3
S_3m(S_3m1==0 & S_3m2==1 & S_3m3==1) = 2; % 2&3
S_3m(biomass3 == 0) = NaN;

%% Interpolate Where Models Don't Use Oxygen

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_o2_e');

% 1 Model
[X,Y] = meshgrid(intlCon1,trsptCon1);
V = cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_1m = double(logical(V)); O_1m(biomass1 == 0) = NaN;

% 2 Models
[X,Y] = meshgrid(intlCon2,trsptCon2);
V = cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_2m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_2m2 = double(logical(V));
O_2m = zeros(numel(trsptCon),numel(intlCon));
O_2m(O_2m1==1 & O_2m2==1) = 2;
O_2m(O_2m1==1 & O_2m2==0) = 1;
O_2m(O_2m1==0 & O_2m2==1) = 1;
O_2m(biomass2 == 0) = NaN;

% 3 Models
[X,Y] = meshgrid(intlCon3,trsptCon3);
V = cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_3m1 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_3m2 = double(logical(V));
V = cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)');
V = interp2(X,Y,V, intlCon_mat,trsptCon_mat, 'linear',0); V(V > 0) = 0;
O_3m3 = double(logical(V));
O_3m = zeros(numel(trsptCon),numel(intlCon));
O_3m(O_3m1==1 & O_3m2==1 & O_3m3==1) = 3; % all
O_3m(O_3m1==1 & O_3m2==0 & O_3m3==0) = 1; % 1 only
O_3m(O_3m1==0 & O_3m2==1 & O_3m3==0) = 1; % 2 only
O_3m(O_3m1==0 & O_3m2==0 & O_3m3==1) = 1; % 3 only
O_3m(O_3m1==1 & O_3m2==1 & O_3m3==0) = 2; % 1&2
O_3m(O_3m1==1 & O_3m2==0 & O_3m3==1) = 2; % 1&3
O_3m(O_3m1==0 & O_3m2==1 & O_3m3==1) = 2; % 2&3
O_3m(biomass3 == 0) = NaN;

%% Parse Metabolite Data

[exchMets_idx,~] = identifyExchMets(Ecoli,model1{1}{1}.exch_idx);
metNames = Ecoli.metNames(exchMets_idx);
% remove medium metabolites
[~,rm_idx,~] = intersect(model1{1}{1}.mets(exchMets_idx),model1{1}{1}.mediumMets);
exchMetsNum_1m(rm_idx) = [];
exchMetsNum_2m(rm_idx) = [];
exchMetsNum_3m(rm_idx) = [];
exchMets_1m(rm_idx) = [];
exchMets_2m(rm_idx) = [];
exchMets_3m(rm_idx) = [];
metNames(rm_idx) = [];
% remove metabolites that are never secreted or exchanged
rm_idx = find(exchMetsNum_1m==0 & exchMetsNum_2m==0 & exchMetsNum_3m==0);
exchMetsNum_1m(rm_idx) = [];
exchMetsNum_2m(rm_idx) = [];
exchMetsNum_3m(rm_idx) = [];
exchMets_1m(rm_idx) = [];
exchMets_2m(rm_idx) = [];
exchMets_3m(rm_idx) = [];
metNames(rm_idx) = [];

clear rm_idx X Y V exchMets_idx

%% Pathways

% Pathways
intl_idx = model2{1}{1}.intl_idx;
pathwayNames = Ecoli.subSystems(intl_idx);
pathwayNames = unique(pathwayNames); 
rm_idx = find(arrayfun(@(x) isempty(pathwayNames{x}),1:numel(pathwayNames))); pathwayNames(rm_idx) = []; clear rm_idx % remove empty
[~,rm_idx,~] = intersect(pathwayNames,'Putative'); pathwayNames(rm_idx) = []; clear rm_idx % remove Putative
[~,rm_idx,~] = intersect(pathwayNames,'Unassigned'); pathwayNames(rm_idx) = []; clear rm_idx % remove Unassigned
[~,rm_idx,~] = intersect(pathwayNames,'Biomass and maintenance functions'); pathwayNames(rm_idx) = []; clear rm_idx % remove Biomass and Maintenance Functions
pathwayNames = strrep(strrep(pathwayNames, 'metabolism','Metabolism'), 'reactions','Reactions');

% Distance Between Models in a Pathway
[X,Y] = meshgrid(intlCon2,trsptCon2);
pathwayJaccDist = cell(numel(pathwayNames),1);
pathwayEucDist = cell(numel(pathwayNames),1);
for pathNum = 1:numel(pathwayNames)
    intlRxns_idx = find(ismember(Ecoli.subSystems,pathwayNames{pathNum}));
    V1 = zeros(numel(trsptCon2),numel(intlCon2));
    V2 = zeros(numel(trsptCon2),numel(intlCon2));
    for trsptNum = 1:numel(trsptCon2)
        for intlNum = find(model2{trsptNum}{1}.biomass)
            [~,intl_idx,~] = intersect(intlCon2,model2{trsptNum}{1}.intl_con(intlNum));
            % Jaccard Distance
            int_2m1 = model2{trsptNum}{1}.int(intlRxns_idx,intlNum);
            int_2m2 = model2{trsptNum}{2}.int(intlRxns_idx,intlNum);
            V1(trsptNum,intl_idx) = pdist2(int_2m1',int_2m2','jaccard');
            % Euclidean Distance
            flux_2m1 = model2{trsptNum}{1}.flux(intlRxns_idx,intlNum);
            flux_2m2 = model2{trsptNum}{2}.flux(intlRxns_idx,intlNum);
            V2(trsptNum,intl_idx) = pdist2(flux_2m1',flux_2m2','euclidean');
        end
    end
    pathwayJaccDist{pathNum} = interp2(X,Y,V1, intlCon_mat,trsptCon_mat, 'linear',0); clear V1
    pathwayJaccDist{pathNum}(isnan(pathwayJaccDist{pathNum})) = 0;
    pathwayJaccDist{pathNum}(biomass2 == 0) = NaN;
    pathwayEucDist{pathNum} = interp2(X,Y,V2, intlCon_mat,trsptCon_mat, 'linear',0); clear V2
    pathwayEucDist{pathNum}(isnan(pathwayEucDist{pathNum})) = 0;
    pathwayEucDist{pathNum}(biomass2 == 0) = NaN;
end

%% Save Data

save(['DOLMN_Parsed/' loadDataName '_interp.mat'], ...
    'x1','x2','x3','y1','y2','y3','z1','z2','z3', ... % biomass flux growth boundary
    'trsptCon','intlCon','biomass*','jaccDist_*','eucDist_*', ... % biomass flux, Jaccard distance, Euclidean distance
    'metNames','exchMets*','numExchMets_*', ... % exchanged metabolites
    'C_*','N_*','P_*','S_*','O_*', ... % where models don't used ...
    'pathwayNames','pathwayJaccDist','pathwayEucDist'); % pathways
