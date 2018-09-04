clear variables
close all
clc

%% Load Data

loadDataName = 'Ecoli_iJR904';

% Load Metabolic Model
load('Ecoli_iJR904_newS.mat')

% Load DOLMN Data
load(['DOLMN_Parsed\iJR904\' loadDataName '.mat'])

%% Biomass

biomass1 = cell2mat(arrayfun(@(x) model1{x}{1}.biomass,1:numel(model1), 'Uni',false)');
biomass2 = cell2mat(arrayfun(@(x) model2{x}{1}.biomass,1:numel(model2), 'Uni',false)');
biomass3 = cell2mat(arrayfun(@(x) model3{x}{1}.biomass,1:numel(model3), 'Uni',false)');

%% Biomass Flux Growth Boundary

% Matrices
biomass1_nan = biomass1; biomass1_nan(biomass1_nan == 0) = NaN;
biomass2_nan = biomass2; biomass2_nan(biomass2_nan == 0) = NaN;
biomass3_nan = biomass3; biomass3_nan(biomass3_nan == 0) = NaN;

% 1 Model
B = bwboundaries(~isnan(biomass1_nan));
x1 = trsptCon1(B{1}(:,1)); y1 = intlCon1(B{1}(:,2));
rm_idx = find(x1 == max(x1) & y1 > min(y1)); x1(rm_idx) = []; y1(rm_idx) = [];
rm_idx = find(y1 == max(y1) & x1 > min(x1)); x1(rm_idx) = []; y1(rm_idx) = [];
idx = find(x1 == max(x1)); x1 = [x1(idx:end), x1(1:idx-1)]; y1 = [y1(idx:end), y1(1:idx-1)];
z1 = arrayfun(@(k) biomass1(trsptCon1 == x1(k),intlCon1 == y1(k)),1:numel(x1));
% 2 Models
B = bwboundaries(~isnan(biomass2_nan));
x2 = trsptCon2(B{1}(:,1)); y2 = intlCon2(B{1}(:,2));
rm_idx = find(x2 == max(x2) & y2 > min(y2)); x2(rm_idx) = []; y2(rm_idx) = [];
rm_idx = find(y2 == max(y2) & x2 > min(x2)); x2(rm_idx) = []; y2(rm_idx) = [];
idx = find(x2 == max(x2)); x2 = [x2(idx:end), x2(1:idx-1)]; y2 = [y2(idx:end), y2(1:idx-1)];
z2 = arrayfun(@(k) biomass2(trsptCon2 == x2(k),intlCon2 == y2(k)),1:numel(x2));
% 3 Models
B = bwboundaries(~isnan(biomass3_nan));
x3 = trsptCon3(B{1}(:,1)); y3 = intlCon3(B{1}(:,2));
rm_idx = find(x3 == max(x3) & y3 > min(y3)); x3(rm_idx) = []; y3(rm_idx) = [];
rm_idx = find(y3 == max(y3) & x3 > min(x3)); x3(rm_idx) = []; y3(rm_idx) = [];
idx = find(x3 == max(x3)); x3 = [x3(idx:end), x3(1:idx-1)]; y3 = [y3(idx:end), y3(1:idx-1)];
z3 = arrayfun(@(k) biomass3(trsptCon3 == x3(k),intlCon3 == y3(k)),1:numel(x3));

x1(2:end) = x1(2:end) - 0.5; x2(2:end) = x2(2:end) - 0.5; x3(2:end) = x3(2:end) - 0.5;
y1(1:end-1) = y1(1:end-1) - 0.5; y2(1:end-1) = y2(1:end-1) - 0.5; y3(1:end-1) = y3(1:end-1) - 0.5;

%% Exchanged/Secreted Metabolites

% 1 Model
exchMetsNum_1m = zeros(Ecoli.N_u,1);
exchMets_1m = cell(Ecoli.N_u,1);
numExchMets_1m = cell2mat(arrayfun(@(x) sum(model1{x}{1}.flux(model1{x}{1}.exch_idx,:) > 0),1:numel(model1), 'Uni',false)');
numExchMets_1m(biomass1 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    exchMets_1m{metNum} = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(model1{x}{1}.exch_idx(metNum),:),1:numel(model1), 'Uni',false)')));
    exchMetsNum_1m(metNum) = sum(exchMets_1m{metNum}(:));
    exchMets_1m{metNum}(exchMets_1m{metNum} == 0) = NaN;
    exchMets_1m{metNum}(biomass1 == 0) = NaN;
end

% 2 Models
exchMetsNum_2m = zeros(Ecoli.N_u,1);
exchMets_2m = cell(Ecoli.N_u,1);
numExchMets_2m = cell2mat(arrayfun(@(x) sum(model2{x}{1}.exchMets.model1_to_model2 + model2{x}{1}.exchMets.model2_to_model1),1:numel(model2), 'Uni',false)');
numExchMets_2m(biomass2 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    exchMets_2m{metNum} = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.exchMets.model1_to_model2(model2{x}{1}.exch_idx(metNum),:) + model2{x}{1}.exchMets.model2_to_model1(model2{x}{1}.exch_idx(metNum),:),1:numel(model2), 'Uni',false)')));
    exchMetsNum_2m(metNum) = sum(exchMets_2m{metNum}(:));
    exchMets_2m{metNum}(exchMets_2m{metNum} == 0) = NaN;
    exchMets_2m{metNum}(biomass2 == 0) = NaN;
end

% 3 Models
exchMetsNum_3m = zeros(Ecoli.N_u,1);
exchMets_3m = cell(Ecoli.N_u,1);
numExchMets_3m = cell2mat(arrayfun(@(x) sum(logical(model3{x}{1}.exchMets.model1_to_model2 + model3{x}{1}.exchMets.model2_to_model1 + ...
    model3{x}{1}.exchMets.model1_to_model3 + model3{x}{1}.exchMets.model3_to_model1 + ...
    model3{x}{1}.exchMets.model2_to_model3 + model3{x}{1}.exchMets.model3_to_model2)),1:numel(model3), 'Uni',false)');
numExchMets_3m(biomass3 == 0) = NaN;
for metNum = 1:Ecoli.N_u
    exchMets_3m{metNum} = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.exchMets.model1_to_model2(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model2_to_model1(model3{x}{1}.exch_idx(metNum),:) + ...
        model3{x}{1}.exchMets.model1_to_model3(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model3_to_model1(model3{x}{1}.exch_idx(metNum),:) + ...
        model3{x}{1}.exchMets.model2_to_model3(model3{x}{1}.exch_idx(metNum),:) + model3{x}{1}.exchMets.model3_to_model2(model3{x}{1}.exch_idx(metNum),:),1:numel(model3), 'Uni',false)')));
    exchMetsNum_3m(metNum) = sum(exchMets_3m{metNum}(:));
    exchMets_3m{metNum}(exchMets_3m{metNum} == 0) = NaN;
    exchMets_3m{metNum}(biomass3 == 0) = NaN;
end

%% Jaccard & Euclidean Distance

% 2 Models - Jaccard
intl_idx = model2{1}{1}.intl_idx;
jaccDist_2m = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.int(intl_idx,:)',model2{x}{2}.int(intl_idx,:)','jaccard')),1:numel(model2), 'Uni',false))';
jaccDist_2m(biomass2 == 0) = NaN;
% 2 Models - Euclidean
eucDist_2m = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.flux(intl_idx,:)',model2{x}{2}.flux(intl_idx,:)','euclidean')),1:numel(model2), 'Uni',false))';
eucDist_2m(biomass2 == 0) = NaN;

% 3 Models - Jaccard
jaccDist_3m12 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.int(intl_idx,:)'),double(model3{x}{2}.int(intl_idx,:)'),'jaccard')),1:numel(model3), 'Uni',false))';
jaccDist_3m13 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.int(intl_idx,:)'),double(model3{x}{3}.int(intl_idx,:)'),'jaccard')),1:numel(model3), 'Uni',false))';
jaccDist_3m23 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{2}.int(intl_idx,:)'),double(model3{x}{3}.int(intl_idx,:)'),'jaccard')),1:numel(model3), 'Uni',false))';
% average jaccard's distance
jaccDist_3m = (jaccDist_3m12 + jaccDist_3m13 + jaccDist_3m23)./3; jaccDist_3m(biomass3 == 0) = NaN;
jaccDist_3m12(biomass3 == 0) = NaN; jaccDist_3m13(biomass3 == 0) = NaN; jaccDist_3m23(biomass3 == 0) = NaN;
% 3 Models - Euclidean
eucDist_3m12 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.flux(intl_idx,:)'),double(model3{x}{2}.flux(intl_idx,:)'),'euclidean')),1:numel(model3), 'Uni',false))';
eucDist_3m13 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{1}.flux(intl_idx,:)'),double(model3{x}{3}.flux(intl_idx,:)'),'euclidean')),1:numel(model3), 'Uni',false))';
eucDist_3m23 = cell2mat(arrayfun(@(x) diag(pdist2(double(model3{x}{2}.flux(intl_idx,:)'),double(model3{x}{3}.flux(intl_idx,:)'),'euclidean')),1:numel(model3), 'Uni',false))';
% average Euclidean distance
eucDist_3m = (eucDist_3m12 + eucDist_3m13 + eucDist_3m23)./3; eucDist_3m(biomass3 == 0) = NaN;
eucDist_3m12(biomass3 == 0) = NaN; eucDist_3m13(biomass3 == 0) = NaN; eucDist_3m23(biomass3 == 0) = NaN;

%% Where Models Don't Use Glucose

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_glc__D_e');

% 1 Model
C_1m = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)')));
C_1m(biomass1 == 0) = NaN;

% 2 Models
C_2m1 = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)')));
C_2m2 = double(logical(cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)')));
C_2m = zeros(numel(trsptCon2),numel(intlCon2));
C_2m(C_2m1==1 & C_2m2==1) = 2;
C_2m(C_2m1==1 & C_2m2==0) = 1;
C_2m(C_2m1==0 & C_2m2==1) = 1;
C_2m(biomass2 == 0) = NaN;

% 3 Models
C_3m1 = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)')));
C_3m2 = double(logical(cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)')));
C_3m3 = double(logical(cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)')));
C_3m = zeros(numel(trsptCon3),numel(intlCon3));
C_3m(C_3m1==1 & C_3m2==1 & C_3m3==1) = 3; % all
C_3m(C_3m1==1 & C_3m2==0 & C_3m3==0) = 1; % 1 only
C_3m(C_3m1==0 & C_3m2==1 & C_3m3==0) = 1; % 2 only
C_3m(C_3m1==0 & C_3m2==0 & C_3m3==1) = 1; % 3 only
C_3m(C_3m1==1 & C_3m2==1 & C_3m3==0) = 2; % 1&2
C_3m(C_3m1==1 & C_3m2==0 & C_3m3==1) = 2; % 1&3
C_3m(C_3m1==0 & C_3m2==1 & C_3m3==1) = 2; % 2&3
C_3m(biomass3 == 0) = NaN;

%% Where Models Don't Use Ammonia

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_nh4_e');

% 1 Model
N_1m = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)')));
N_1m(biomass1 == 0) = NaN;

% 2 Models
N_2m1 = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)')));
N_2m2 = double(logical(cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)')));
N_2m = zeros(numel(trsptCon2),numel(intlCon2));
N_2m(N_2m1==1 & N_2m2==1) = 2;
N_2m(N_2m1==1 & N_2m2==0) = 1;
N_2m(N_2m1==0 & N_2m2==1) = 1;
N_2m(biomass2 == 0) = NaN;

% 3 Models
N_3m1 = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)')));
N_3m2 = double(logical(cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)')));
N_3m3 = double(logical(cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)')));
N_3m = zeros(numel(trsptCon3),numel(intlCon3));
N_3m(N_3m1==1 & N_3m2==1 & N_3m3==1) = 3; % all
N_3m(N_3m1==1 & N_3m2==0 & N_3m3==0) = 1; % 1 only
N_3m(N_3m1==0 & N_3m2==1 & N_3m3==0) = 1; % 2 only
N_3m(N_3m1==0 & N_3m2==0 & N_3m3==1) = 1; % 3 only
N_3m(N_3m1==1 & N_3m2==1 & N_3m3==0) = 2; % 1&2
N_3m(N_3m1==1 & N_3m2==0 & N_3m3==1) = 2; % 1&3
N_3m(N_3m1==0 & N_3m2==1 & N_3m3==1) = 2; % 2&3
N_3m(biomass3 == 0) = NaN;

%% Where Models Don't Use Phosphate

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_pi_e');

% 1 Model
P_1m = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)')));
P_1m(biomass1 == 0) = NaN;

% 2 Models
P_2m1 = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)')));
P_2m2 = double(logical(cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)')));
P_2m = zeros(numel(trsptCon2),numel(intlCon2));
P_2m(P_2m1==1 & P_2m2==1) = 2;
P_2m(P_2m1==1 & P_2m2==0) = 1;
P_2m(P_2m1==0 & P_2m2==1) = 1;
P_2m(biomass2 == 0) = NaN;

% 3 Models
P_3m1 = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)')));
P_3m2 = double(logical(cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)')));
P_3m3 = double(logical(cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)')));
P_3m = zeros(numel(trsptCon3),numel(intlCon3));
P_3m(P_3m1==1 & P_3m2==1 & P_3m3==1) = 3; % all
P_3m(P_3m1==1 & P_3m2==0 & P_3m3==0) = 1; % 1 only
P_3m(P_3m1==0 & P_3m2==1 & P_3m3==0) = 1; % 2 only
P_3m(P_3m1==0 & P_3m2==0 & P_3m3==1) = 1; % 3 only
P_3m(P_3m1==1 & P_3m2==1 & P_3m3==0) = 2; % 1&2
P_3m(P_3m1==1 & P_3m2==0 & P_3m3==1) = 2; % 1&3
P_3m(P_3m1==0 & P_3m2==1 & P_3m3==1) = 2; % 2&3
P_3m(biomass3 == 0) = NaN;

%% Where Models Don't Use Sulfate

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_so4_e');

% 1 Model
S_1m = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)')));
S_1m(biomass1 == 0) = NaN;

% 2 Models
S_2m1 = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)')));
S_2m2 = double(logical(cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)')));
S_2m = zeros(numel(trsptCon2),numel(intlCon2));
S_2m(S_2m1==1 & S_2m2==1) = 2;
S_2m(S_2m1==1 & S_2m2==0) = 1;
S_2m(S_2m1==0 & S_2m2==1) = 1;
S_2m(biomass2 == 0) = NaN;

% 3 Models
S_3m1 = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)')));
S_3m2 = double(logical(cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)')));
S_3m3 = double(logical(cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)')));
S_3m = zeros(numel(trsptCon3),numel(intlCon3));
S_3m(S_3m1==1 & S_3m2==1 & S_3m3==1) = 3; % all
S_3m(S_3m1==1 & S_3m2==0 & S_3m3==0) = 1; % 1 only
S_3m(S_3m1==0 & S_3m2==1 & S_3m3==0) = 1; % 2 only
S_3m(S_3m1==0 & S_3m2==0 & S_3m3==1) = 1; % 3 only
S_3m(S_3m1==1 & S_3m2==1 & S_3m3==0) = 2; % 1&2
S_3m(S_3m1==1 & S_3m2==0 & S_3m3==1) = 2; % 1&3
S_3m(S_3m1==0 & S_3m2==1 & S_3m3==1) = 2; % 2&3
S_3m(biomass3 == 0) = NaN;

%% Where Models Don't Use Oxygen

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_o2_e');

% 1 Model
O_1m = double(logical(cell2mat(arrayfun(@(x) model1{x}{1}.flux(idx,:),1:numel(model1), 'Uni',false)')));
O_1m(biomass1 == 0) = NaN;

% 2 Models
O_2m1 = double(logical(cell2mat(arrayfun(@(x) model2{x}{1}.flux(idx,:),1:numel(model2), 'Uni',false)')));
O_2m2 = double(logical(cell2mat(arrayfun(@(x) model2{x}{2}.flux(idx,:),1:numel(model2), 'Uni',false)')));
O_2m = zeros(numel(trsptCon2),numel(intlCon2));
O_2m(O_2m1==1 & O_2m2==1) = 2;
O_2m(O_2m1==1 & O_2m2==0) = 1;
O_2m(O_2m1==0 & O_2m2==1) = 1;
O_2m(biomass2 == 0) = NaN;

% 3 Models
O_3m1 = double(logical(cell2mat(arrayfun(@(x) model3{x}{1}.flux(idx,:),1:numel(model3), 'Uni',false)')));
O_3m2 = double(logical(cell2mat(arrayfun(@(x) model3{x}{2}.flux(idx,:),1:numel(model3), 'Uni',false)')));
O_3m3 = double(logical(cell2mat(arrayfun(@(x) model3{x}{3}.flux(idx,:),1:numel(model3), 'Uni',false)')));
O_3m = zeros(numel(trsptCon3),numel(intlCon3));
O_3m(O_3m1==1 & O_3m2==1 & O_3m3==1) = 3; % all
O_3m(O_3m1==1 & O_3m2==0 & O_3m3==0) = 1; % 1 only
O_3m(O_3m1==0 & O_3m2==1 & O_3m3==0) = 1; % 2 only
O_3m(O_3m1==0 & O_3m2==0 & O_3m3==1) = 1; % 3 only
O_3m(O_3m1==1 & O_3m2==1 & O_3m3==0) = 2; % 1&2
O_3m(O_3m1==1 & O_3m2==0 & O_3m3==1) = 2; % 1&3
O_3m(O_3m1==0 & O_3m2==1 & O_3m3==1) = 2; % 2&3
O_3m(biomass3 == 0) = NaN;

%% Exchange Flux Correlations
% Rho: Spearman's Rho, measures the strength of association between continuous variables (ranking)
%   +: model* secretes/takes up both metabolites
%   -: model* secretes one metabolite and takes up the other metabolite
% Phi: Special case of Pearson's correlation coefficient, measures the degree of association between BINARY variables
%   +: metabolite A is exchanged at the same constraints as metabolite B
%   -: metabolite A is exchanged at the constraints that metabolite B is not exchanged

[exchMets_idx,~] = identifyExchMets(Ecoli,model1{1}{1}.exch_idx);
exchMetNames = Ecoli.metNames(exchMets_idx);
[~,idx,~] = intersect(model1{1}{1}.mets(exchMets_idx),model1{1}{1}.mediumMets);
mediumMets = model1{1}{1}.metNames(exchMets_idx(idx)); clear idx
data_flux_1m  = cell2mat(arrayfun(@(x) model1{x}{1}.flux(model1{x}{1}.exch_idx(exchMets_idx),:), 1:numel(model1), 'Uni',false));
data_flux_2m1 = cell2mat(arrayfun(@(x) model2{x}{1}.flux(model2{x}{1}.exch_idx(exchMets_idx),:), 1:numel(model2), 'Uni',false));
data_flux_2m2 = cell2mat(arrayfun(@(x) model2{x}{2}.flux(model2{x}{2}.exch_idx(exchMets_idx),:), 1:numel(model2), 'Uni',false));
data_flux_3m1 = cell2mat(arrayfun(@(x) model3{x}{1}.flux(model3{x}{1}.exch_idx(exchMets_idx),:), 1:numel(model3), 'Uni',false));
data_flux_3m2 = cell2mat(arrayfun(@(x) model3{x}{2}.flux(model3{x}{2}.exch_idx(exchMets_idx),:), 1:numel(model3), 'Uni',false));
data_flux_3m3 = cell2mat(arrayfun(@(x) model3{x}{3}.flux(model3{x}{3}.exch_idx(exchMets_idx),:), 1:numel(model3), 'Uni',false));

% Remove No-Growth Cases
rm_idx1 = find(biomass1(:) == 0);
data_flux_1m(:,rm_idx1) = [];
rm_idx2 = find(biomass2(:) == 0);
data_flux_2m1(:,rm_idx2) = []; data_flux_2m2(:,rm_idx2) = [];
rm_idx3 = find(biomass3(:) == 0);
data_flux_3m1(:,rm_idx3) = []; data_flux_3m2(:,rm_idx3) = []; data_flux_3m3(:,rm_idx3) = [];

% 1 Model
[rho1,~] = corr(data_flux_1m', 'Type','Spearman'); rho1(isnan(rho1)) = 0;
Z = linkage(rho1,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P1_rho] = dendrogram(Z,0); close gcf
temp = cell2mat(arrayfun(@(x) exchMets_1m{x}(:),1:numel(exchMets_1m), 'Uni',false)); temp(isnan(temp)) = 0; temp(rm_idx1,:) = [];
[phi1,~] = corr(temp, 'Type','Pearson'); phi1(isnan(phi1)) = 0; clear temp
Z = linkage(phi1,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P1_phi] = dendrogram(Z,0); close gcf

% 2 Models
[rho2,~] = corr([data_flux_2m1, data_flux_2m2]', 'Type','Spearman'); rho2(isnan(rho2)) = 0;
Z = linkage(rho2,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P2_rho] = dendrogram(Z,0); close gcf
temp = cell2mat(arrayfun(@(x) exchMets_2m{x}(:),1:numel(exchMets_2m), 'Uni',false)); temp(isnan(temp)) = 0; temp(rm_idx2,:) = [];
[phi2,~] = corr(temp, 'Type','Pearson'); phi2(isnan(phi2)) = 0; clear temp
Z = linkage(phi2,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P2_phi] = dendrogram(Z,0); close gcf

% 3 Models
[rho3,~] = corr([data_flux_3m1, data_flux_3m2, data_flux_3m3]', 'Type','Spearman'); rho3(isnan(rho3)) = 0;
Z = linkage(rho3,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P3_rho] = dendrogram(Z,0); close gcf
temp = cell2mat(arrayfun(@(x) exchMets_3m{x}(:),1:numel(exchMets_3m), 'Uni',false)); temp(isnan(temp)) = 0; temp(rm_idx3,:) = [];
[phi3,~] = corr(temp, 'Type','Pearson'); phi3(isnan(phi3)) = 0; clear temp
Z = linkage(phi3,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P3_phi] = dendrogram(Z,0); close gcf

%% Parse Metabolite Data

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

%% PCA

intl_idx = model2{1}{1}.intl_idx;
rxns = model2{1}{1}.rxns(intl_idx);

% Where 1 Model Stops Growing
temp = arrayfun(@(x) min(model1{x}{1}.intl_con(find(model1{x}{1}.biomass)))-1,1:numel(model1));
min_intlCon1 = interp1(trsptCon1,temp,trsptCon2); clear temp

% 1 Model
flux1m = cell2mat(arrayfun(@(x) model1{x}{1}.flux(intl_idx,:),1:numel(model1), 'Uni',false));
Nt1m = cell2mat(arrayfun(@(x) model1{x}{1}.trspt_con,1:numel(model1), 'Uni',false));
Ni1m = cell2mat(arrayfun(@(x) model1{x}{1}.intl_con,1:numel(model1), 'Uni',false));
% 2 Models
flux2m1 = cell2mat(arrayfun(@(x) model2{x}{1}.flux(intl_idx,:),1:numel(model2), 'Uni',false));
flux2m2 = cell2mat(arrayfun(@(x) model2{x}{2}.flux(intl_idx,:),1:numel(model2), 'Uni',false));
Nt2m = cell2mat(arrayfun(@(x) model2{x}{1}.trspt_con,1:numel(model2), 'Uni',false));
Ni2m = cell2mat(arrayfun(@(x) model2{x}{1}.intl_con,1:numel(model2), 'Uni',false));

% Remove All Zero Cases
% 1 Model
rm_idx = find(ismember(flux1m',zeros(1,size(flux1m,1)),'rows'));
flux1m(:,rm_idx) = []; Nt1m(rm_idx) = []; Ni1m(rm_idx) = []; clear rm_idx
% 2 Models
rm_idx2m1 = find(ismember(flux2m1',zeros(1,size(flux2m1,1)),'rows'));
rm_idx2m2 = find(ismember(flux2m2',zeros(1,size(flux2m1,1)),'rows'));
rm_idx = intersect(rm_idx2m1,rm_idx2m2);
flux2m1(:,rm_idx) = []; flux2m2(:,rm_idx) = []; Nt2m(rm_idx) = []; Ni2m(rm_idx) = []; clear rm_idx

% Biomass Flux
[~,bio_idx,~] = intersect(rxns,'BIOMASS_Ecoli'); rxns(bio_idx) = [];
bioFlux1m  =  flux1m(bio_idx,:);  flux1m(bio_idx,:) = [];
bioFlux2m1 = flux2m1(bio_idx,:); flux2m1(bio_idx,:) = [];
bioFlux2m2 = flux2m2(bio_idx,:); flux2m2(bio_idx,:) = [];

% PCA
flux = [flux1m./repmat(bioFlux1m,numel(rxns),1), ...
    flux2m1./repmat(bioFlux2m1,numel(rxns),1), ...
    flux2m2./repmat(bioFlux2m2,numel(rxns),1)]; % normalize by biomass flux
[coeff,score,~,~,explained] = pca(flux');
score1m  = score(1:size(flux1m,2),:);
score2m1 = score(size(flux1m,2)+1:size(flux1m,2)+size(flux2m1,2),:);
score2m2 = score(size(flux1m,2)+size(flux2m1,2)+1:end,:);
clear score

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
% 2 Models
pathwayJaccDist = cell(numel(pathwayNames),1);
pathwayEucDist = cell(numel(pathwayNames),1);
for pathNum = 1:numel(pathwayNames)
    intlRxns_idx = find(ismember(Ecoli.subSystems,pathwayNames{pathNum}));
    pathwayJaccDist{pathNum} = zeros(size(biomass2));
    pathwayEucDist{pathNum} = zeros(size(biomass2));
    for trsptNum = 1:numel(trsptCon2)-1
        for intlNum = find(model2{trsptNum}{1}.biomass)
            [~,intl_idx,~] = intersect(intlCon2,model2{trsptNum}{1}.intl_con(intlNum));
            % Jaccard Distance
            int_2m1 = model2{trsptNum}{1}.int(intlRxns_idx,intlNum);
            int_2m2 = model2{trsptNum}{2}.int(intlRxns_idx,intlNum);
            pathwayJaccDist{pathNum}(trsptNum,intl_idx) = pdist2(int_2m1',int_2m2','jaccard');
            % Euclidean Distance
            flux_2m1 = model2{trsptNum}{1}.flux(intlRxns_idx,intlNum);
            flux_2m2 = model2{trsptNum}{2}.flux(intlRxns_idx,intlNum);
            pathwayEucDist{pathNum}(trsptNum,intl_idx) = pdist2(flux_2m1',flux_2m2','euclidean');
        end
    end
    pathwayJaccDist{pathNum}(isnan(pathwayJaccDist{pathNum})) = 0;
    pathwayJaccDist{pathNum}(biomass2 == 0) = NaN;
    pathwayEucDist{pathNum}(isnan(pathwayEucDist{pathNum})) = 0;
    pathwayEucDist{pathNum}(biomass2 == 0) = NaN;
end

%% Pathways and Exchanged Metabolites
% Point-Biserial Correlation: Special case of Pearson's correlation coefficient, measures the degree of association between 
%       a continuous variable and a binary variable
%   +: metabolite A is exchanged at the same constraints as pathway B is distributed
%   -: metabolite A is exchanged at the constraints that pathways B is not distributed

% Correlation
rm_idx = find(biomass2(:) == 0); % remove no-growth cases
corr_pathMet_jacc = zeros(numel(pathwayNames),numel(metNames));
corr_pathMet_euc = zeros(numel(pathwayNames),numel(metNames));
for pathNum = 1:numel(pathwayNames)
    for metNum = 1:numel(metNames)
        X = pathwayJaccDist{pathNum}(:); X(isnan(X)) = 0;
        Y = exchMets_2m{metNum}(:);  Y(isnan(Y)) = 0;
        X(rm_idx) = []; Y(rm_idx) = [];
        corr_pathMet_jacc(pathNum,metNum) = corr(X,Y, 'Type','Pearson', 'Rows','Complete'); clear X Y
        
        X = pathwayEucDist{pathNum}(:); X(isnan(X)) = 0;
        Y = exchMets_2m{metNum}(:);  Y(isnan(Y)) = 0;
        X(rm_idx) = []; Y(rm_idx) = [];
        corr_pathMet_euc(pathNum,metNum) = corr(X,Y, 'Type','Pearson', 'Rows','Complete'); clear X Y
    end
end
corr_pathMet_jacc(isnan(corr_pathMet_jacc)) = 0; clear rm_idx
corr_pathMet_euc(isnan(corr_pathMet_euc)) = 0; clear rm_idx

% Hierarchical Clustering
Z = linkage(corr_pathMet_jacc,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P_pathJacc] = dendrogram(Z,0); close gcf
Z = linkage(corr_pathMet_jacc','single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P_metsJacc] = dendrogram(Z,0); close gcf

% Hierarchical Clustering
Z = linkage(corr_pathMet_euc,'single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P_pathEuc] = dendrogram(Z,0); close gcf
Z = linkage(corr_pathMet_euc','single','euclidean'); % get tree of hierarchical clusters (cluster closest distances)
figure; [~, ~, P_metsEuc] = dendrogram(Z,0); close gcf

%% Reactions Kept

intlRxns = model1{1}{1}.rxns(model1{1}{1}.intl_idx);
trsptRxns = model1{1}{1}.rxns(model1{1}{1}.trspt_idx);

% 1 Model
minNi = min(intlCon1(arrayfun(@(x) find(biomass1(x,:) == 0,1,'last'),1:size(biomass1,1))))+1;
min_intlRxns1m = []; min_trsptRxns1m = [];
intlRxns1m = []; trsptRxns1m = [];
% Minimum Transport Reactions
for intlNum = find(model1{1}{1}.biomass)
    min_trsptRxns1m = [min_trsptRxns1m, trsptRxns(logical(model1{1}{1}.int(model1{1}{1}.trspt_idx,intlNum)))];
end
for trsptNum = 1:numel(model1)
    for intlNum = find(model1{trsptNum}{1}.biomass)
        % Minimum Intracellular Reactions
        if model1{trsptNum}{1}.intl_con(intlNum) == minNi
            min_intlRxns1m = [min_intlRxns1m, intlRxns(logical(model1{trsptNum}{1}.int(model1{trsptNum}{1}.intl_idx,intlNum)))];
        end
        % Intracellular Reactions
        intlRxns1m = [intlRxns1m; intlRxns(logical(model1{trsptNum}{1}.int(model1{trsptNum}{1}.intl_idx,intlNum)))];
        % Transport Reactions
        trsptRxns1m = [trsptRxns1m; trsptRxns(logical(model1{trsptNum}{1}.int(model1{trsptNum}{1}.trspt_idx,intlNum)))];
    end
end

% 2 Models
minNi = min(intlCon2(arrayfun(@(x) find(biomass2(x,:) == 0,1,'last'),1:size(biomass2,1))))+1;
min_intlRxns2m1 = []; min_intlRxns2m2 = [];
min_trsptRxns2m1 = []; min_trsptRxns2m2 = [];
intlRxns2m1 = []; trsptRxns2m1 = [];
intlRxns2m2 = []; trsptRxns2m2 = [];
% Minimum Transport Reactions
for intlNum = find(model2{1}{1}.biomass)
    min_trsptRxns2m1 = [min_trsptRxns2m1, trsptRxns(logical(model2{1}{1}.int(model2{1}{1}.trspt_idx,intlNum)))];
    min_trsptRxns2m2 = [min_trsptRxns2m2, trsptRxns(logical(model2{1}{2}.int(model2{1}{2}.trspt_idx,intlNum)))];
end
for trsptNum = 1:numel(model2)
    for intlNum = find(model2{trsptNum}{1}.biomass)
        % Minimum Intracellular Reactions
        if model2{trsptNum}{1}.intl_con(intlNum) == minNi
            min_intlRxns2m1 = [min_intlRxns2m1, intlRxns(logical(model2{trsptNum}{1}.int(model2{trsptNum}{1}.intl_idx,intlNum)))];
            min_intlRxns2m2 = [min_intlRxns2m2, intlRxns(logical(model2{trsptNum}{2}.int(model2{trsptNum}{2}.intl_idx,intlNum)))];
        end
        % Intracellular Reactions
        intlRxns2m1 = [intlRxns2m1; intlRxns(logical(model2{trsptNum}{1}.int(model2{trsptNum}{1}.intl_idx,intlNum)))];
        intlRxns2m2 = [intlRxns2m2; intlRxns(logical(model2{trsptNum}{2}.int(model2{trsptNum}{2}.intl_idx,intlNum)))];
        % Transport Reactions
        trsptRxns2m1 = [trsptRxns2m1; trsptRxns(logical(model2{trsptNum}{1}.int(model2{trsptNum}{1}.trspt_idx,intlNum)))];
        trsptRxns2m2 = [trsptRxns2m2; trsptRxns(logical(model2{trsptNum}{2}.int(model2{trsptNum}{2}.trspt_idx,intlNum)))];
    end
end

% 3 Models
minNi = min(intlCon3(arrayfun(@(x) find(biomass3(x,:) == 0,1,'last'),1:size(biomass3,1))))+1;
min_intlRxns3m1 = []; min_intlRxns3m2 = []; min_intlRxns3m3 = [];
min_trsptRxns3m1 = []; min_trsptRxns3m2 = []; min_trsptRxns3m3 = [];
intlRxns3m1 = []; trsptRxns3m1 = [];
intlRxns3m2 = []; trsptRxns3m2 = [];
intlRxns3m3 = []; trsptRxns3m3 = [];
% Minimum Transport Reactions
for intlNum = find(model3{1}{1}.biomass)
    min_trsptRxns3m1 = [min_trsptRxns3m1, trsptRxns(logical(model3{1}{1}.int(model3{1}{1}.trspt_idx,intlNum)))];
    min_trsptRxns3m2 = [min_trsptRxns3m2, trsptRxns(logical(model3{1}{2}.int(model3{1}{2}.trspt_idx,intlNum)))];
    min_trsptRxns3m3 = [min_trsptRxns3m3, trsptRxns(logical(model3{1}{3}.int(model3{1}{2}.trspt_idx,intlNum)))];
end
for trsptNum = 1:numel(model3)
    for intlNum = find(model3{trsptNum}{1}.biomass)
        % Minimum Intracellular Reactions
        if model3{trsptNum}{1}.intl_con(intlNum) == minNi
            min_intlRxns3m1 = [min_intlRxns3m1, intlRxns(logical(model3{trsptNum}{1}.int(model3{trsptNum}{1}.intl_idx,intlNum)))];
            min_intlRxns3m2 = [min_intlRxns3m2, intlRxns(logical(model3{trsptNum}{2}.int(model3{trsptNum}{2}.intl_idx,intlNum)))];
            min_intlRxns3m3 = [min_intlRxns3m3, intlRxns(logical(model3{trsptNum}{3}.int(model3{trsptNum}{3}.intl_idx,intlNum)))];
        end
        % Intracellular Reactions
        intlRxns3m1 = [intlRxns3m1; intlRxns(logical(model3{trsptNum}{1}.int(model3{trsptNum}{1}.intl_idx,intlNum)))];
        intlRxns3m2 = [intlRxns3m2; intlRxns(logical(model3{trsptNum}{2}.int(model3{trsptNum}{2}.intl_idx,intlNum)))];
        intlRxns3m3 = [intlRxns3m3; intlRxns(logical(model3{trsptNum}{3}.int(model3{trsptNum}{3}.intl_idx,intlNum)))];
        % Transport Reactions
        trsptRxns3m1 = [trsptRxns3m1; trsptRxns(logical(model3{trsptNum}{1}.int(model3{trsptNum}{1}.trspt_idx,intlNum)))];
        trsptRxns3m2 = [trsptRxns3m2; trsptRxns(logical(model3{trsptNum}{2}.int(model3{trsptNum}{2}.trspt_idx,intlNum)))];
        trsptRxns3m3 = [trsptRxns3m3; trsptRxns(logical(model3{trsptNum}{3}.int(model3{trsptNum}{3}.trspt_idx,intlNum)))];
    end
end

%% Min Intracellular Reactions
clc

% 1 Model
min_intlRxns1m_same = arrayfun(@(x) numel(min_intlRxns1m(:,x)),1:size(min_intlRxns1m,2));
C = arrayfun(@(x) min_intlRxns1m(:,x),1:size(min_intlRxns1m,2),'Uni',false);
intlRxns1m_same = C{1};
for ii = 2:numel(C)
    intlRxns1m_same = intersect(intlRxns1m_same,C{ii});
end
% Make the Same Size
for ii = 1:numel(C)
    n = max(min_intlRxns1m_same) - numel(C{ii});
    N = numel(C{ii});
    for jj = 1:n
        C{ii}(N+jj) = {''};
    end
end
intlRxns1m_diff = setxor(unique([C{:}]),intlRxns1m_same);
intlRxns1m_diff(ismember(intlRxns1m_diff,'')) = [];

% 2 Models
min_intlRxns2m_same = arrayfun(@(x) numel(intersect(min_intlRxns2m1(:,x),min_intlRxns2m2(:,x))),1:size(min_intlRxns2m1,2));
C = arrayfun(@(x) intersect(min_intlRxns2m1(:,x),min_intlRxns2m2(:,x)),1:size(min_intlRxns2m1,2),'Uni',false);
intlRxns2m_same = C{1};
for ii = 2:numel(C)
    intlRxns2m_same = intersect(intlRxns2m_same,C{ii});
end
% Make the Same Size
for ii = 1:numel(C)
    n = max(min_intlRxns2m_same) - numel(C{ii});
    N = numel(C{ii});
    for jj = 1:n
        C{ii}(N+jj) = {''};
    end
end
intlRxns2m_diff = setxor(unique([C{:}]),intlRxns2m_same);
intlRxns2m_diff(ismember(intlRxns2m_diff,'')) = [];

% 3 Models
min_intlRxns3m_same = arrayfun(@(x) numel(intersect(intersect(min_intlRxns3m1(:,x),min_intlRxns3m2(:,x)),min_intlRxns3m3(:,x))),1:size(min_intlRxns2m1,2));
C = arrayfun(@(x) intersect(intersect(min_intlRxns3m1(:,x),min_intlRxns3m2(:,x)),min_intlRxns3m3(:,x)),1:size(min_intlRxns3m1,2),'Uni',false);
intlRxns3m_same = C{1};
for ii = 2:numel(C)
    intlRxns3m_same = intersect(intlRxns3m_same,C{ii});
end
% Make the Same Size
for ii = 1:numel(C)
    n = max(min_intlRxns3m_same) - numel(C{ii});
    N = numel(C{ii});
    for jj = 1:n
        C{ii}(N+jj) = {''};
    end
end
intlRxns3m_diff = setxor(unique([C{:}]),intlRxns3m_same);
intlRxns3m_diff(ismember(intlRxns3m_diff,'')) = [];

% Min Transport Reactions
min_trsptRxns = union(union(unique(min_trsptRxns1m), ...
    union(unique(min_trsptRxns2m1),unique(min_trsptRxns2m2))), ...
    union(union(unique(min_trsptRxns3m1),unique(min_trsptRxns3m2)),unique(min_trsptRxns3m3)));
min_trsptRxns1m = histcounts(categorical(min_trsptRxns1m),categorical(min_trsptRxns));
min_trsptRxns2m1 = histcounts(categorical(min_trsptRxns2m1),categorical(min_trsptRxns));
min_trsptRxns2m2 = histcounts(categorical(min_trsptRxns2m2),categorical(min_trsptRxns));
min_trsptRxns3m1 = histcounts(categorical(min_trsptRxns3m1),categorical(min_trsptRxns));
min_trsptRxns3m2 = histcounts(categorical(min_trsptRxns3m2),categorical(min_trsptRxns));
min_trsptRxns3m3 = histcounts(categorical(min_trsptRxns3m3),categorical(min_trsptRxns));

%% Save Data

save(['DOLMN_Parsed/' loadDataName '_summary.mat'], ...
    'x1','x2','x3','y1','y2','y3','z1','z2','z3', ... % biomass flux growth boundary
    'trsptCon*','intlCon*','biomass*','jaccDist_*','eucDist_*', ... % biomass flux, Jaccard distance, Euclidean distance
    'metNames','exchMets*','numExchMets_*','corr_pathMet*', ... % exchanged metabolites
    'rho*','P*','phi*','exchMets_idx','exchMetNames','mediumMets', ... % exchanged metabolites
    'C_*','N_*','P_*','S_*','O_*', ... % where models don't used ...
    'coeff','score*','explained','Nt1m','Ni1m','Nt2m','Ni2m', ... % PCA
    'pathwayNames','pathwayJaccDist','pathwayEucDist', ... % pathways
    'min_intlRxns*','min_trsptRxns*','intlRxns*','trsptRxns*'); % reactions kept

