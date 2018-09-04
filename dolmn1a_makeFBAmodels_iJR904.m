clear variables
close all
clc

%% Load Data

if ~isdir('DOLMN_Parsed/iJR904/'); mkdir('DOLMN_Parsed/iJR904/'); end

saveDataName = 'Ecoli_iJR904';
model_title = 'E. coli iJR904';

% Load Metabolic Model
load('Ecoli_iJR904_newS.mat')

% Medium
[~,exch_rxns] = identifyExchRxns(Ecoli);
[metInds,~] = identifyMetsFromRxns(Ecoli,exch_rxns);
medium = Ecoli.mets(metInds);
medium_rxns = Ecoli.rxns(exch_rxns);

% ATPM Reaction
atpm_name = 'ATPM';

% Transport Reaction Constraints
trsptCon1 = [9:22, 24, 34, 38, 205];
trsptCon2_3 = [9:45, 205];

% Biomass Flux Lower Bound
bio_lb = 0.1;

%% Make FBA Models

% Each Transport Constraint in 1, 2, & 3 Models
trsptCon = intersect(trsptCon1,trsptCon2_3);
for trsptCon_num = trsptCon
    disp(['T_trspt=' int2str(trsptCon_num)])
    
    % 1 Model
    disp('     1 Model')
    data1 = load(['DOLMN_Output\iJR904\CVX_iJR904_K1_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data1.trspt_con = trsptCon_num.*ones(size(data1.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data1.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data1.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data1.biomass(:,nan_idx) = 0;
        data1.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data1.model)
            data1.model{model_num}.biomass(nan_idx) = 0;
            data1.model{model_num}.flux(:,nan_idx) = 0;
            data1.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data1.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    
    % Biomass Flux Lower Bound
    data1.bio_lb = bio_lb;
    % Process
    [models1,base_model1,model_flux1] = algorithm2models(data1,Ecoli,[model_title ' 1 Model'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_1Model_TrsptCon_' int2str(trsptCon_num) '.mat'],'models1','base_model1','model_flux1')
    
    % 2 Models
    disp('     2 Models')
    data2 = load(['DOLMN_Output\iJR904\CVX_iJR904_K2_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data2.trspt_con = trsptCon_num.*ones(size(data2.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data2.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data2.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data2.biomass(:,nan_idx) = 0;
        data2.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data2.model)
            data2.model{model_num}.biomass(nan_idx) = 0;
            data2.model{model_num}.flux(:,nan_idx) = 0;
            data2.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data2.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    % Biomass Flux Lower Bound
    data2.bio_lb = bio_lb;
    % Process
    [models2,base_model2,model_flux2] = algorithm2models(data2,Ecoli,[model_title ' 2 Models'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_2Models_TrsptCon_' int2str(trsptCon_num) '.mat'],'models2','base_model2','model_flux2')
    
    % 3 Models
    disp('     3 Models')
    data3 = load(['DOLMN_Output\iJR904\CVX_iJR904_K3_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data3.trspt_con = trsptCon_num.*ones(size(data3.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data3.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data3.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data3.biomass(:,nan_idx) = 0;
        data3.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data3.model)
            data3.model{model_num}.biomass(nan_idx) = 0;
            data3.model{model_num}.flux(:,nan_idx) = 0;
            data3.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data3.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    % Biomass Flux Lower Bound
    data3.bio_lb = bio_lb;
    % Process
    [models3,base_model3,model_flux3] = algorithm2models(data3,Ecoli,[model_title ' 3 Models'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_3Models_TrsptCon_' int2str(trsptCon_num) '.mat'],'models3','base_model3','model_flux3')
    
    clear data1* data2* data3* models* base_model* model_flux* opt_status*
    close all
    disp('=========================')
end

% Each Transport Constraint in 1 Model
trsptCon1only = setdiff(trsptCon1,trsptCon);
for trsptCon_num = trsptCon1only
    disp(['T_trspt=' int2str(trsptCon_num)])
    
    % 1 Model
    disp('     1 Model')
    data1 = load(['DOLMN_Output\iJR904\CVX_iJR904_K1_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data1.trspt_con = trsptCon_num.*ones(size(data1.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data1.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data1.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data1.biomass(:,nan_idx) = 0;
        data1.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data1.model)
            data1.model{model_num}.biomass(nan_idx) = 0;
            data1.model{model_num}.flux(:,nan_idx) = 0;
            data1.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data1.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    % Biomass Flux Lower Bound
    data1.bio_lb = bio_lb;
    % Process
    [models1,base_model1,model_flux1] = algorithm2models(data1,Ecoli,[model_title ' 1 Model'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_1Model_TrsptCon_' int2str(trsptCon_num) '.mat'],'models1','base_model1','model_flux1')
    
    clear data1* models* base_model* model_flux* opt_status*
    close all
    disp('=========================')
end

% Each Transport Constraint in 2 & 3 Models
trsptCon2_3only = setdiff(trsptCon2_3,trsptCon);
for trsptCon_num = trsptCon2_3only
    disp(['T_trspt=' int2str(trsptCon_num)])

    % 2 Models
    disp('     2 Models')
    data2 = load(['DOLMN_Output\iJR904\CVX_iJR904_K2_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data2.trspt_con = trsptCon_num.*ones(size(data2.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data2.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data2.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data2.biomass(:,nan_idx) = 0;
        data2.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data2.model)
            data2.model{model_num}.biomass(nan_idx) = 0;
            data2.model{model_num}.flux(:,nan_idx) = 0;
            data2.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data2.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    % Biomass Flux Lower Bound
    data2.bio_lb = bio_lb;
    % Process
    [models2,base_model2,model_flux2] = algorithm2models(data2,Ecoli,[model_title ' 2 Models'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_2Models_TrsptCon_' int2str(trsptCon_num) '.mat'],'models2','base_model2','model_flux2')
    
    % 3 Models
    disp('     3 Models')
    data3 = load(['DOLMN_Output\iJR904\CVX_iJR904_K3_Nt' int2str(trsptCon_num) '.mat']);
    % Data
    data3.trspt_con = trsptCon_num.*ones(size(data3.sparse_con));
    % Check for Infeasible Solutions
    [~,nan_idx] = find(isnan(data3.biomass)); nan_idx = unique(nan_idx);
    [~,bio_idx] = find(data3.biomass==0); bio_idx = min(bio_idx);
    if nan_idx == bio_idx-1
        data3.biomass(:,nan_idx) = 0;
        data3.flux_norm(:,nan_idx) = 0;
        for model_num = 1:numel(data3.model)
            data3.model{model_num}.biomass(nan_idx) = 0;
            data3.model{model_num}.flux(:,nan_idx) = 0;
            data3.model{model_num}.int(:,nan_idx) = 0;
        end
    elseif ~isempty(nan_idx)
        error('Error in data: Infeasible solution is not at boundary')
    end
    % Check for Numerical Problems
    [~,num_idx] = find(data3.biomass > 1); num_idx = unique(num_idx);
    if ~isempty(num_idx)
        error('Error in data: Numerical problem')
    end
    % Biomass Flux Lower Bound
    data3.bio_lb = bio_lb;
    % Process
    [models3,base_model3,model_flux3] = algorithm2models(data3,Ecoli,[model_title ' 3 Models'],medium_rxns,atpm_name);
    % Save Data
    save(['DOLMN_Parsed\iJR904\' saveDataName '_3Models_TrsptCon_' int2str(trsptCon_num) '.mat'],'models3','base_model3','model_flux3')
    
    clear data2* data3* models* base_model* model_flux* opt_status*
    close all
    disp('=========================')
end

%% Combine Data

% 1 Model
model1 = cell(size(trsptCon1)); % Pre-Allocate
for trsptCon_num = 1:numel(trsptCon1)
    % Load Data
    load(['DOLMN_Parsed\iJR904\' saveDataName '_1Model_TrsptCon_' int2str(trsptCon1(trsptCon_num)) '.mat'],'model_flux1')
    % Save to New Struct
    model1{trsptCon_num} = model_flux1; clear model_flux1
end

% 2 Models
model2 = cell(size(trsptCon2_3)); % Pre-Allocate
for trsptCon_num = 1:numel(trsptCon2_3)
    % Load Data
    load(['DOLMN_Parsed\iJR904\' saveDataName '_2Models_TrsptCon_' int2str(trsptCon2_3(trsptCon_num)) '.mat'],'model_flux2')
    % Save to New Struct
    model2{trsptCon_num} = model_flux2; clear model_flux2
end

% 3 Models
model3 = cell(size(trsptCon2_3)); % Pre-Allocate
for trsptCon_num = 1:numel(trsptCon2_3)
    % Load Data
    load(['DOLMN_Parsed\iJR904\' saveDataName '_3Models_TrsptCon_' int2str(trsptCon2_3(trsptCon_num)) '.mat'],'model_flux3')
    % Save to New Struct
    model3{trsptCon_num} = model_flux3; clear model_flux3
end

%% Pad because added a column whenever biomass didn't reach zero

% 1 Model
num_intlCon1 = arrayfun(@(x) numel(model1{x}{1}.intl_con), 1:numel(model1));
for ii = find(num_intlCon1 < max(num_intlCon1))
    model1{ii}{1}.intl_con = [model1{ii}{1}.intl_con(1)-1, model1{ii}{1}.intl_con];
    model1{ii}{1}.trspt_con = [model1{ii}{1}.trspt_con(1), model1{ii}{1}.trspt_con];
    model1{ii}{1}.biomass = [0, model1{ii}{1}.biomass];
    model1{ii}{1}.flux = [zeros(numel(model1{ii}{1}.rxns),1), model1{ii}{1}.flux];
    model1{ii}{1}.int = [zeros(numel(model1{ii}{1}.rxns),1), model1{ii}{1}.int];
    model1{ii}{1}.totalExchFlux = [zeros(numel(model1{ii}{1}.exch_idx),1), model1{ii}{1}.totalExchFlux];
    model1{ii}{1}.warningFlag.BinaryVariableFlux = [0, model1{ii}{1}.warningFlag.BinaryVariableFlux];
    model1{ii}{1}.warningFlag.CalcExchFlux = [0, model1{ii}{1}.warningFlag.CalcExchFlux];
    model1{ii}{1}.warningFlag.FBA_sol = [0, model1{ii}{1}.warningFlag.FBA_sol];
    model1{ii}{1}.warningFlag.sol = [0, model1{ii}{1}.warningFlag.sol];
    model1{ii}{1}.warningFlag.Sv_viol = [0, model1{ii}{1}.warningFlag.Sv_viol];
    model1{ii}{1}.warningFlag.Sv = [zeros(size(model1{1}{1}.mets)), model1{ii}{1}.warningFlag.Sv];
    model1{ii}{1}.warningFlag.bounds_viol = [0, model1{ii}{1}.warningFlag.bounds_viol];
    model1{ii}{1}.warningFlag.intl_bounds = [zeros(size(model1{1}{1}.intl_idx)), model1{ii}{1}.warningFlag.intl_bounds];
end

% 2 Models
num_intlCon2 = arrayfun(@(x) numel(model2{x}{1}.intl_con), 1:numel(model2));
for ii = find(num_intlCon2 < max(num_intlCon2))
    for model_num = 1:2
        model2{ii}{model_num}.intl_con = [model2{ii}{model_num}.intl_con(1)-1, model2{ii}{model_num}.intl_con];
        model2{ii}{model_num}.trspt_con = [model2{ii}{model_num}.trspt_con(1), model2{ii}{model_num}.trspt_con];
        model2{ii}{model_num}.biomass = [0, model2{ii}{model_num}.biomass];
        model2{ii}{model_num}.flux = [zeros(numel(model2{ii}{model_num}.rxns),1), model2{ii}{model_num}.flux];
        model2{ii}{model_num}.int = [zeros(numel(model2{ii}{model_num}.rxns),1), model2{ii}{model_num}.int];
        model2{ii}{model_num}.totalExchFlux = [zeros(numel(model2{ii}{model_num}.exch_idx),1), model2{ii}{model_num}.totalExchFlux];
        model2{ii}{model_num}.exchMets.fluxType12 = [zeros(numel(model2{ii}{model_num}.exch_idx),1), model2{ii}{model_num}.exchMets.fluxType12];
        model2{ii}{model_num}.exchMets.model1_to_model2 = [zeros(numel(model2{ii}{model_num}.exch_idx),1), model2{ii}{model_num}.exchMets.model1_to_model2];
        model2{ii}{model_num}.exchMets.model2_to_model1 = [zeros(numel(model2{ii}{model_num}.exch_idx),1), model2{ii}{model_num}.exchMets.model2_to_model1];
        model2{ii}{model_num}.warningFlag.BinaryVariableFlux = [0, model2{ii}{model_num}.warningFlag.BinaryVariableFlux];
        model2{ii}{model_num}.warningFlag.CalcExchFlux = [0, model2{ii}{model_num}.warningFlag.CalcExchFlux];
        model2{ii}{model_num}.warningFlag.FBA_sol = [0, model2{ii}{model_num}.warningFlag.FBA_sol];
        model2{ii}{model_num}.warningFlag.sol = [0, model2{ii}{model_num}.warningFlag.sol];
        model2{ii}{model_num}.warningFlag.Sv_viol = [0, model2{ii}{model_num}.warningFlag.Sv_viol];
        model2{ii}{model_num}.warningFlag.Sv = [zeros(size(model2{1}{1}.mets)), model2{ii}{model_num}.warningFlag.Sv];
        model2{ii}{model_num}.warningFlag.bounds_viol = [0, model2{ii}{model_num}.warningFlag.bounds_viol];
        model2{ii}{model_num}.warningFlag.intl_bounds = [zeros(size(model2{1}{1}.intl_idx)), model2{ii}{model_num}.warningFlag.intl_bounds];
    end
end

% 3 Models
num_intlCon3 = arrayfun(@(x) numel(model3{x}{1}.intl_con), 1:numel(model3));
for ii = find(num_intlCon3 < max(num_intlCon3))
    for model_num = 1:3
        model3{ii}{model_num}.intl_con = [model3{ii}{model_num}.intl_con(1)-1, model3{ii}{model_num}.intl_con];
        model3{ii}{model_num}.trspt_con = [model3{ii}{model_num}.trspt_con(1), model3{ii}{model_num}.trspt_con];
        model3{ii}{model_num}.biomass = [0, model3{ii}{model_num}.biomass];
        model3{ii}{model_num}.flux = [zeros(numel(model3{ii}{model_num}.rxns),1), model3{ii}{model_num}.flux];
        model3{ii}{model_num}.int = [zeros(numel(model3{ii}{model_num}.rxns),1), model3{ii}{model_num}.int];
        model3{ii}{model_num}.totalExchFlux = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.totalExchFlux];
        model3{ii}{model_num}.exchMets.fluxType12 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.fluxType12];
        model3{ii}{model_num}.exchMets.model1_to_model2 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model1_to_model2];
        model3{ii}{model_num}.exchMets.model2_to_model1 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model2_to_model1];
        model3{ii}{model_num}.exchMets.fluxType13 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.fluxType13];
        model3{ii}{model_num}.exchMets.model1_to_model3 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model1_to_model3];
        model3{ii}{model_num}.exchMets.model3_to_model1 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model3_to_model1];
        model3{ii}{model_num}.exchMets.fluxType23 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.fluxType23];
        model3{ii}{model_num}.exchMets.model2_to_model3 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model2_to_model3];
        model3{ii}{model_num}.exchMets.model3_to_model2 = [zeros(numel(model3{ii}{model_num}.exch_idx),1), model3{ii}{model_num}.exchMets.model3_to_model2];
        model3{ii}{model_num}.warningFlag.BinaryVariableFlux = [0, model3{ii}{model_num}.warningFlag.BinaryVariableFlux];
        model3{ii}{model_num}.warningFlag.CalcExchFlux = [0, model3{ii}{model_num}.warningFlag.CalcExchFlux];
        model3{ii}{model_num}.warningFlag.FBA_sol = [0, model3{ii}{model_num}.warningFlag.FBA_sol];
        model3{ii}{model_num}.warningFlag.sol = [0, model3{ii}{model_num}.warningFlag.sol];
        model3{ii}{model_num}.warningFlag.Sv_viol = [0, model3{ii}{model_num}.warningFlag.Sv_viol];
        model3{ii}{model_num}.warningFlag.Sv = [zeros(size(model3{1}{1}.mets)), model3{ii}{model_num}.warningFlag.Sv];
        model3{ii}{model_num}.warningFlag.bounds_viol = [0, model3{ii}{model_num}.warningFlag.bounds_viol];
        model3{ii}{model_num}.warningFlag.intl_bounds = [zeros(size(model3{1}{1}.intl_idx)), model3{ii}{model_num}.warningFlag.intl_bounds];
    end
end

%% Save Data

% Sparsity Constraints
intlCon1 = model1{end}{1}.intl_con; % intracellular sparsity constraints
trsptCon1 = arrayfun(@(x) model1{x}{1}.trspt_con(1), 1:numel(model1)); % transport sparsity constraints
intlCon2 = model2{end}{1}.intl_con; % intracellular sparsity constraints
trsptCon2 = arrayfun(@(x) model2{x}{1}.trspt_con(1), 1:numel(model2)); % transport sparsity constraints
intlCon3 = model3{end}{1}.intl_con; % intracellular sparsity constraints
trsptCon3 = arrayfun(@(x) model3{x}{1}.trspt_con(1), 1:numel(model3)); % transport sparsity constraints

% Save
save(['DOLMN_Parsed\iJR904\' saveDataName '.mat'],...
    'model1','model2','model3','trsptCon1','trsptCon2','trsptCon3','intlCon1','intlCon2','intlCon3')


