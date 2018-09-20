clear variables
close all
clc

%% Load Data

loadDataName = 'Ecoli_Core';

% Load Metabolic Model
load('Ecoli_core_newS.mat')

% Load DOLMN Data
load(fullfile('DOLMN_Parsed','Core',[loadDataName '.mat']))

%% Biomass

biomass1 = cell2mat(arrayfun(@(x) model1{x}{1}.biomass,1:numel(model1), 'Uni',false)');
biomass2 = cell2mat(arrayfun(@(x) model2{x}{1}.biomass,1:numel(model2), 'Uni',false)');

%% Biomass Flux Growth Boundary

% Matrices
biomass1_nan = biomass1; biomass1_nan(biomass1_nan == 0) = NaN;
biomass2_nan = biomass2; biomass2_nan(biomass2_nan == 0) = NaN;

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

x1(2:end) = x1(2:end) - 0.5; x2(2:end) = x2(2:end) - 0.5;
y1(1:end-1) = y1(1:end-1) - 0.5; y2(1:end-1) = y2(1:end-1) - 0.5;

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

%% Jaccard & Euclidean Distance

% 2 Models - Jaccard
intl_idx = model2{1}{1}.intl_idx;
jaccDist_2m = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.int(intl_idx,:)',model2{x}{2}.int(intl_idx,:)','jaccard')),1:numel(model2), 'Uni',false))';
jaccDist_2m(biomass2 == 0) = NaN;
% 2 Models - Euclidean
eucDist_2m = cell2mat(arrayfun(@(x) diag(pdist2(model2{x}{1}.flux(intl_idx,:)',model2{x}{2}.flux(intl_idx,:)','euclidean')),1:numel(model2), 'Uni',false))';
eucDist_2m(biomass2 == 0) = NaN;

%% Where Models Don't Use Glucose

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_glc(e)');

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

%% Where Models Don't Use Ammonia

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_nh4(e)');

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

%% Where Models Don't Use Phosphate

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_pi(e)');

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

%% Where Models Don't Use Oxygen

[~,idx,~] = intersect(model1{1}{1}.rxns,'EX_o2(e)');

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

% Remove No-Growth Cases
rm_idx1 = find(biomass1(:) == 0);
data_flux_1m(:,rm_idx1) = [];
rm_idx2 = find(biomass2(:) == 0);
data_flux_2m1(:,rm_idx2) = []; data_flux_2m2(:,rm_idx2) = [];

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

%% Parse Metabolite Data

metNames = Ecoli.metNames(exchMets_idx);
% remove medium metabolites
[~,rm_idx,~] = intersect(model1{1}{1}.mets(exchMets_idx),model1{1}{1}.mediumMets);
exchMetsNum_1m(rm_idx) = [];
exchMetsNum_2m(rm_idx) = [];
exchMets_1m(rm_idx) = [];
exchMets_2m(rm_idx) = [];
metNames(rm_idx) = [];
% remove metabolites that are never secreted or exchanged
rm_idx = find(exchMetsNum_1m==0 & exchMetsNum_2m==0);
exchMetsNum_1m(rm_idx) = [];
exchMetsNum_2m(rm_idx) = [];
exchMets_1m(rm_idx) = [];
exchMets_2m(rm_idx) = [];
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
[~,bio_idx,~] = intersect(rxns,'Biomass_Ecoli_core_w_GAM'); rxns(bio_idx) = [];
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

%% Split TCA

tcaRxns = {'ACONTa';'ACONTb';'AKGDH';'CS';'FRD7';'FUM';'ICDHyr';'MDH';'SUCDi';'SUCOAS'}; % 'MALS'; 'ICL'
tcaJaccDist = zeros(size(biomass2));
tcaEucDist = zeros(size(biomass2));
for trsptNum = 1:numel(trsptCon2)-1
    [~,tcaRxns_idx,~] = intersect(Ecoli.rxns,tcaRxns);
    for intlNum = find(model2{trsptNum}{1}.biomass)
        [~,intl_idx,~] = intersect(intlCon2,model2{trsptNum}{1}.intl_con(intlNum));
        % Jaccard Distance
        int_2m1 = double(model2{trsptNum}{1}.int(tcaRxns_idx,intlNum));
        int_2m2 = double(model2{trsptNum}{2}.int(tcaRxns_idx,intlNum));
        tcaJaccDist(trsptNum,intl_idx) = pdist2(int_2m1',int_2m2','jaccard');
        % Euclidean Distance
        flux_2m1 = model2{trsptNum}{1}.flux(tcaRxns_idx,intlNum);
        flux_2m2 = model2{trsptNum}{2}.flux(tcaRxns_idx,intlNum);
        tcaEucDist(trsptNum,intl_idx) = pdist2(flux_2m1',flux_2m2','euclidean');
    end
end
tcaJaccDist(biomass2 == 0) = NaN;
tcaEucDist(biomass2 == 0) = NaN;

%% Save Data

save(fullfile('DOLMN_Parsed',[loadDataName '_summary.mat']), ...
    'x1','x2','y1','y2','z1','z2', ... % biomass flux growth boundary
    'trsptCon*','intlCon*','biomass*','jaccDist_*','eucDist_*', ... % biomass flux, Jaccard distance, Euclidean distance
    'metNames','exchMets*','numExchMets_*', ... % exchanged metabolites
    'rho*','P*','phi*','exchMets_idx','exchMetNames','mediumMets', ... % exchanged metabolites
    'C_*','N_*','P_*','O_*', ... % where models don't use ...
    'coeff','score*','explained','Nt1m','Ni1m','Nt2m','Ni2m', ... % PCA
    'tcaJaccDist','tcaEucDist'); % TCA

