%% input folder core_input:
%   S matrix, lb,ub,metNames,rxnNames and essential reactions
%% output folder core_output_1/core_output_2:
%   1 or 2 core biomass decrease for stricter internal Constraint
%   	if Trspt Constraint is fixed number in file names.
%%
% this code is faster since the optimal solutions are provided as start
% point. But without good start, it may take hours.

% require gurobi and cvx
%   http://cvxr.com/cvx/doc/gurobi.html
% in original code, transport reaction is called exchange reaction so
% 	variable names may be confusing
% 20170524 t has same length as x
% 1,2 core equal_biomass, biomass>0.1 for each bacteria

function []=GITecoli_core_equal12(N_ecoli,save_output1,save_output2,N_cpu)

load('core_input/ecoli_core_newS.mat','newS','mets','rxns', 'rxnNames','N_u','N_e','N_i','M_i','M_e','lb','ub','c');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust bound to speed up
x_lb=lb;
x_lb(c==1)=0.1; %biomass>0.1 for each bacteria
unique_x_lb=unique(x_lb)'; % -100.0000  -10.0000         0    8.3900
id_ATPM=find(ismember(x_lb,x_lb(x_lb> x_lb(c==1) ) ) );
% rxns( id_ATPM )% 'ATPM'

x_ub=ub;
Max = 100;
x_lb(x_lb<-Max)=-Max;
x_ub(x_ub>Max)=Max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sx=0 ==> S(:,-53)x(-53)=-S(:,53)x(53) due to x(53)==8.39
%name also need to change!!
% from 1 S to community S
NUM=N_ecoli;%2%3;% number of bacteria
newS_store=newS;% used for b_e

newS(:,id_ATPM)=[];
x_lb(id_ATPM)=[]; x_ub(id_ATPM)=[];
c(id_ATPM)=[];
rxns(id_ATPM)=[];
N_i=N_i-1;

Su=sparse( M_e+NUM*M_i,N_u+NUM*(N_i+N_e) );
Su(1:M_e, 1:N_u)=newS(1:M_e, 1:N_u);
newS_e1=newS(1:M_e, 1+N_u:end);
Su(1:M_e, 1+N_u:end)=repmat(newS_e1,1,NUM);
newS_e2_i=newS(1+M_e:end, 1+N_u:end);
Su(1+M_e:end, 1+N_u:end)=kron(eye(NUM),newS_e2_i);

[m,n] = size(Su);
x_lb2=[x_lb(1:N_u); repmat(x_lb(N_u+1:end),NUM,1)];
x_ub2=[x_ub(1:N_u);repmat(x_ub(N_u+1:end),NUM,1)];
c=[c(1:N_u);repmat(c(N_u+1:end),NUM,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = zeros((2*n),1);
obj(1:length(c))=c;
biomass_id=find(obj);
%equality constraints
A_e = sparse(m,2*n);
A_e(1:m,1:n) = Su;
temp = full( -newS_store(:,id_ATPM)*lb(id_ATPM));% one ATMP
b_e =[temp(1:M_e); repmat(temp(M_e+1:end),NUM,1)];

% equal_biomass constraints
if NUM>1
    A_e2 = sparse(NUM-1,n+n);
    b_e2 =zeros(NUM-1,1);
    for k=2:NUM
        A_e2(k-1,biomass_id(1))=1;
        A_e2(k-1,biomass_id(k))=-1;
    end
    A_e =[A_e ;A_e2 ];
    b_e =[b_e ;b_e2 ];
end

% inequality constraints
A_ineq_1 = sparse(2*n,2*n);
A_ineq_1(:, 1:n) = [eye(n); -eye(n)];
A_ineq_1(:, n+1:2*n) = [-diag(x_ub2); diag(x_lb2)];
b1_ineq = zeros(2*n,1);

% specify the type of variables
Ctype = char('C'*ones(2*n,1))';
Ctype(end - n+1:end) = char('B'*ones(n,1))';

% lower and upper bounds of variables in the problem
LB = [x_lb2;  zeros(n,1)];
UB = [x_ub2;  ones(n,1)];
%LB(n+1:end)=-Inf; UB(n+1:end)=Inf;
LB(n+1:n+N_u)=1;

% essential_rxns
if NUM==1
    load('core_input/core_1_essential_rxns01.mat')%../
    LB(n+N_u+N_e +atleast1)=1;
else
    load('core_input/core_2_essential_rxns01_equal_biomass.mat')
    LB(n+N_u+N_e +atleast2)=1;
    LB(n+N_u+N_e +N_i+N_e +atleast2)=1;
end

% % regulation constraints on T
T_ineq = sparse(NUM+NUM,2*n);
% overall bound t_min <= sum t_i <= t_max
id_con=n+N_u+N_e +1:n+N_u+N_e +N_i;
for k=1:NUM
    T_ineq (k,id_con) =1;
    id_con=id_con+N_e +N_i;
end
%exchange sparsity constraints
id_con=n+N_u +1:n+N_u+N_e;
for k=NUM+1:NUM+NUM
    T_ineq (k,id_con) =1;
    id_con=id_con+N_e +N_i;
end

%%do not leave space for n+N_u+N_e+atleast1_exclude_atleast2(i)
%essential_rxns atleast1_exclude_atleast2
if NUM>1
    T_ineq2 = sparse(length(atleast1_exclude_atleast2),n+n);
    for i=1:length(atleast1_exclude_atleast2)
        T_ineq2(i,[n+N_u+N_e+atleast1_exclude_atleast2(i),n+N_u+N_e+N_i+N_e+atleast1_exclude_atleast2(i)]) =-1;
    end
    t_ineq2=-ones(length(atleast1_exclude_atleast2),1);
end

if NUM==1
    sparse_con=[N_i,32:-1:20];
    sparse_EX=[N_e,12:-1:9];
else
    sparse_con=[N_i,30:-1:15];
    sparse_EX=[N_e,16:-1:9];
end



for  kk=1:length(sparse_EX)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flux_all=zeros(n+n,size(sparse_con,2));
    %     %     if NUM==1
    %     %         filename = ['core_primal1_sparse_EX_20170524_flux_all/core_primal' num2str( NUM ) '_sparse_EX_' num2str( sparse_EX(kk) ) '_bio01_20170524_flux_all.mat']
    %     %     else
    %     %         filename = ['core_primal2_sparse_EX_20170524_flux_all/core_primal' num2str( NUM ) '_sparse_EX_' num2str( sparse_EX(kk)  ) '_equal_bio01_20170524_flux_all.mat']
    %     %     end
    filename = ['core_output_' num2str( NUM ) '/core_' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk)  ) '.mat'];

    flux_allold=load(filename,'flux_all');
    flux_all(:,1:size(flux_allold.flux_all,2))=flux_allold.flux_all;
    flux_norm=zeros(3,length(sparse_con));
    
    if NUM==1
        biomass=zeros(1,length(sparse_con));
        biomassold=flux_all(obj>0,:);
    else
        biomass=zeros(NUM+1,length(sparse_con));
        biomassold=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
    end
    biomass(:,1:length(biomassold))=biomassold;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flux_all_row_id=cell(size(flux_all,1),1);
    flux_all_row_id(1:N_u)={'U'};
    flux_all_row_id(n+1:n+N_u)={'TU'};
    for k=1:NUM
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+1:N_u+(N_e+N_i)*(k-1)+N_e)={['E' num2str( k)]};
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+N_e+1:N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['I' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+1:n+N_u+(N_e+N_i)*(k-1)+N_e)={['TE' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+N_e+1:n+N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['TI' num2str( k)]};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     filename0 = ['core_primal' num2str( NUM ) '_sparse_EX_' num2str( sparse_EX(kk-1) ) '_equal_bio01_20170524_flux_all.mat'];
    %     start=load(filename0,'flux_all');
    fprintf('1st MIP optimization to max biomass \n');
    start=load(filename,'flux_all');
    for k=1:length(sparse_con)
        %         if sum(start.flux_all(obj>0,k))==0
        %             break
        %         end
        % upper and lowber bounds of sum_i t_i
        t_a_max = sparse_con(k);
        t_ineq =[repmat(t_a_max,NUM,1);repmat(sparse_EX(kk),NUM,1)];
        if NUM==1
            A_ineq = [A_ineq_1; T_ineq];
            b_ineq = [b1_ineq;t_ineq];
        else
            A_ineq = [A_ineq_1; T_ineq; T_ineq2];
            b_ineq = [b1_ineq;  t_ineq; t_ineq2];
        end
        if k>5
            A_ineq=[A_ineq; obj'];
            b_ineq =[b_ineq; biomass(end,k-1)];
        end
        
        clear model;
        model.obj = obj;
        model.A = sparse([A_e; A_ineq]);
        model.sense = [ char('='*ones(size(A_e,1),1));char('<'*ones(size(A_ineq,1),1))];
        model.rhs = [b_e;b_ineq];
        model.vtype = Ctype;
        model.lb = LB;
        model.ub = UB;
        model.modelsense = 'max'; % maximize objective function
        %%
        model.start = nan(size(model.A,2),1);
        model.start(n+N_u+1:2*n) =(abs(start.flux_all(N_u+1:n,k))>10^-7);
        if k>1
            model.start(n+N_u+1:n+N_u+N_e) =(abs(flux_all(N_u+1:N_u+N_e,k-1))>10^-7);
        end
        %gurobi_write(model, 'core2.lp');
        %%
        clear params;
        params.outputflag = 0;
        %   params.Presolve = 2;
        params.TimeLimit = 10;
        params.IntFeasTol=1e-9;%Default value:	1e-5
        %params.MIPGapAbs=0;%Default value:	1e-10
        % params.MIPGap=1e-9;
        params.threads =N_cpu;
        params.FeasibilityTol=1e-9;%Default value:	1e-6
        result = gurobi(model,params);
        
        if  isfield(result,'x') && biomass(end,k)<result.objval
            x =result.x;
            flux_all(:,k)=result.x;
            if NUM==1
                biomass=flux_all(obj>0,:);
            else
                biomass=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
            end
            %             sparsity_in_x=[sum(abs(flux_all(ismember(flux_all_row_id,{'I1'}),:))>10^-7);sum(abs(flux_all(ismember(flux_all_row_id,{'I2'}),:))>10^-7)];
            flux_norm(:,k)=[norm(x(1+N_u:end),1);norm(x,1);norm(x(1:N_u),1)];
        else
            break%continue;%
        end
    end
    %%
    clear model;
    for k=1:NUM
        model{k}.sparse_con = sparse_con;
        model{k}.rxns = rxns;%rxns([exch_flux, trspt_flux1, intl_flux1]);
        model{k}.biomass = biomass(k,:);
        model{k}.flux = flux_all(ismember(flux_all_row_id,{'U',['E' num2str( k)],['I' num2str( k)]}),:);%M1_flux;  sum(ismember(flux_all_row_id,{'U','E1','I1'}))
        model{k}.int =  flux_all(ismember(flux_all_row_id,{'TU',['TE' num2str( k)],['TI' num2str( k)]}),:);
    end
    %     sparsity_ex_x=[sum(abs(flux_all(ismember(flux_all_row_id,{'E1'}),:))>10^-9);sum(abs(flux_all(ismember(flux_all_row_id,{'E2'}),:))>10^-9)];
    %     sparsity_in_x=[sum(abs(flux_all(ismember(flux_all_row_id,{'I1'}),:))>10^-7);sum(abs(flux_all(ismember(flux_all_row_id,{'I2'}),:))>10^-7)];
    fprintf('%d core biomass decrease for stricter internalCon as follows if TrsptCon is %d.\n',NUM,sparse_EX(kk));
    biomass(end,:)
    flux_all(n+1:n+N_u,:)=abs(flux_all(1:N_u,:)>10^-8);
    if save_output1==1
        save(filename,'model','sparse_EX','sparse_con','flux_all','rxns','mets','biomass','biomass_id','flux_norm','flux_all_row_id')
    end
    
    %% CVX
    fprintf('2nd optimization to min flux norm \n');
    
    opt=load(filename);
    flux_allold=load(filename,'flux_all');
    flux_all(:,1:size(flux_allold.flux_all,2))=flux_allold.flux_all;
    flux_norm=zeros(3,length(sparse_con));
    for k=1:length(sparse_con)
        if opt.biomass(1,k)>0
            t=opt.flux_all(n+1:end,k);
            t(1:N_u)=1;
            cvx_clear
            %             cvx_solver gurobi
            cvx_precision best %high
            cvx_begin  quiet
            variable x(n)
            %variable t(n) binary
            minimize( norm(x(N_u+1:end),1) )
            subject to
            Su*x == b_e(1:m);
            x_lb2 <= x <= x_ub2;
            c'*x == opt.biomass(1,k)*NUM;
            diag(x_lb2)*t<= x <= diag(x_ub2)*t;
            for kkk=1:NUM
                x(biomass_id(1))==x(biomass_id(kkk));%equal biomass
            end
            cvx_end
            flux_all(:,k)=[x;t];
            if NUM==1
                biomass=flux_all(obj>0,:);
            else
                biomass=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
            end
            flux_norm(:,k)=[norm(x(1+N_u:end),1);norm(x,1);norm(x(1:N_u),1)];
        else
            break
        end
    end
    
    clear model;
    for k=1:NUM
        model{k}.sparse_con = sparse_con;
        model{k}.rxns = rxns;
        model{k}.biomass = biomass(k,:);
        model{k}.flux = flux_all(ismember(flux_all_row_id,{'U',['E' num2str( k)],['I' num2str( k)]}),:);%M1_flux;  sum(ismember(flux_all_row_id,{'U','E1','I1'}))
        model{k}.int =  flux_all(ismember(flux_all_row_id,{'TU',['TE' num2str( k)],['TI' num2str( k)]}),:);
    end
    flux_all(n+1:n+N_u,:)=abs(flux_all(1:N_u,:)>10^-8);
    %%
    biomass(end,:)
    filenameCVX = ['core_output_' num2str( NUM ) '/CVX_core_' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk)  ) '.mat'];
    
    if save_output2==1
        save(filenameCVX,'model','sparse_EX','sparse_con','flux_all','rxns','mets','biomass','biomass_id','flux_norm','flux_all_row_id')
    end
    
    if biomass(1,1)==0
        break
    end
end
end
