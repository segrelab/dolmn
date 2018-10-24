%% input folder iJR904_input:
%   S matrix, lb,ub,metNames,rxnNames and essential reactions
%% output folder iJR904_output_1/iJR904_output_2/iJR904_output_3:
%   1 or 2 or 3 iJR904 biomass decrease for stricter internal Constraint
%   	if Trspt Constraint is fixed number in file names.
%%
% this code is faster since the optimal solutions are provided as start
% point. But without good start, it may take hours.

% require gurobi and cvx
%   http://cvxr.com/cvx/doc/gurobi.html
% in original code, transport reaction is called exchange reaction so
% 	variable names may be confusing
% t has same length as x
%  biomass>0.1 for each bacteria
% clc;
% clear;

% function iJR904_x_t_equal_random(range)
% range_sparseEX=1:37;
% EXACT=-1 for NUM==3 GAME heuristic solution from NUM-1==2 results (offer initial solutions for NUM>1)
% EXACT=0 for heuristic but fast, EXACT=1 for exact solution but slow
function []=GITecoli_iJR904_equal123(N_ecoli,save_output1,save_output2,N_cpu,EXACT)
% EXACT=0
for repeat=1:(1+20*(EXACT==0))
    tic
    fprintf('heuristic optimization times');
    disp(repeat)
    load('iJR904_input/iJR904_newS20171006.mat','newS','mets','rxns','N_u','N_e','N_i','M_i','M_e','lb','ub','c');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adjust bound to speed up
    x_lb=lb;%
    x_lb(c==1)=0.1; % biomass>0.1 for each bacteria
%     unique_x_lb=unique(x_lb)';
    id_ATPM=find(ismember(x_lb,x_lb(x_lb> x_lb(c==1) ) ) );
    %     rxns( id_ATPM );% 'ATPM'
    x_ub=ub;
    Max = 100;
    x_lb(x_lb<-Max)=-Max;
    x_ub(x_ub>Max)=Max;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sx=0 ==> S(:,-53)x(-53)=-S(:,53)x(53) due to x(53)==8.39
    % from 1 S to community S
    NUM=N_ecoli;%3;% number of bacteria
    newS_store=newS;
    
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
    % equality constraints
    A_e = sparse(m,2*n);
    A_e(1:m,1:n) = Su;
    temp = full( -newS_store(:,id_ATPM)*lb(id_ATPM));
    b_e =[temp(1:M_e); repmat(temp(M_e+1:end),NUM,1)];
    %equal_biomass constraints
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
    %% not important for uptake rxns here
    LB(n+1:n+N_u)=1;
    %essential_rxns
    if NUM==1
        load('iJR904_input/iJR904_x_t_equal_atleast1_essential_rxns_01.mat')
        LB(n+atleast1)=1;
    elseif NUM==2
        load('iJR904_input/iJR904_x_t_equal_atleast1_2_essential_rxns_01.mat')
        LB(n+atleast2)=1;
        LB(n+(N_i+N_e)+atleast2)=1;
    else
        load('iJR904_input/iJR904_x_t_equal_atleast1_3_essential_rxns_01.mat')
        LB(n+atleast3)=1;
        LB(n+(N_i+N_e)+atleast3)=1;
        LB(n+(N_i+N_e)*2+atleast3)=1;
    end
    % % regulation constraints on T
    T_ineq = sparse(NUM+NUM,2*n);
    % overall bound t_min <= sum t_i <= t_max
    % intracellular sparsity constraints
    id_con=n+N_u+N_e +1:n+N_u+N_e +N_i;
    for k=1:NUM
        T_ineq (k,id_con) =1;
        id_con=id_con+N_e +N_i;
    end
    % transport sparsity constraints
    id_con=n+N_u +1:n+N_u+N_e;
    for k=NUM+1:NUM+NUM
        T_ineq (k,id_con) =1;
        id_con=id_con+N_e+N_i;
    end
    % essential_rxns atleast1_exclude_atleast2
    if NUM==2
        len=length(atleast1_exclude_atleast2);
        T_ineq2 = sparse(len,n+n);
        for i=1:len
            id_con=n+atleast1_exclude_atleast2(i);
            for j=1:NUM
                T_ineq2(i,id_con) =-1;
                id_con=id_con+N_e+N_i;
            end
        end
        t_ineq2=-ones(len,1);
    end
    if NUM==3
        len=length(atleast1_exclude_atleast3);
        T_ineq2 = sparse(len,n+n);
        for i=1:len
            id_con=n+atleast1_exclude_atleast3(i);
            for j=1:NUM
                T_ineq2(i,id_con) =-1;
                id_con=id_con+N_e +N_i;
            end
        end
        t_ineq2=-ones(len,1);
    end
    %% constraints on intracellular and transport reactions
    if NUM==1
        sparse_con=[N_i,296:-1:250];
        sparse_EX=[N_e,38,34,24,22:-1:8];
    elseif NUM==2
        sparse_con=[N_i,285:-1:214];
        sparse_EX=[N_e,45:-1:8];
    else
        sparse_con=[N_i,285:-1:194];
        sparse_EX=[N_e,45:-1:8];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flux_all_row_id=cell(n+n,1);
    flux_all_row_id(1:N_u)={'U'};
    flux_all_row_id(n+1:n+N_u)={'TU'};
    for k=1:NUM
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+1:N_u+(N_e+N_i)*(k-1)+N_e)={['E' num2str( k)]};
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+N_e+1:N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['I' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+1:n+N_u+(N_e+N_i)*(k-1)+N_e)={['TE' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+N_e+1:n+N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['TI' num2str( k)]};
    end
    %% loop over transport reactions constraints
    for  kk=1:length(sparse_EX)-1
        flux_all=zeros(n+n,size(sparse_con,2));
%         filename = ['iJR904_output_' num2str( NUM ) '/CVX_iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];
        filename = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];

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
        filename0 = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(max(kk-1,1)) ) '.mat'];
        filename2 = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(min(kk+1,38)) ) '.mat'];
        
        next=load(filename2,'flux_all');
        pre=load(filename0,'flux_all');
        %% loop over intracellular reactions constrains
        for   k=1:length(sparse_con)
%             if biomass(end,max(k-1,1))==0 %|| biomass(end,k)==0
%                 break
%             end
%             %% FOCUS ON POSSIBLE SOLUTION
%             if    biomass(end,k)>0 
%                 continue
%             end
%             if   biomass(end,max(k-1,1))<0.4
%                 continue
%             end
            %%
            if    biomass(end,k)>max(biomass(end,1:k))-0.0001%0.0001 %|| (biomass(end,k)-biomass0.biomass(end,k))>-0.01%(k>1 && biomass(end,k)==biomass(end,k-1))
                continue
            end
            if    biomass(end,k)>biomass(end,max(k-1,1))-0.01%0.0001 %|| (biomass(end,k)-biomass0.biomass(end,k))>-0.01%(k>1 && biomass(end,k)==biomass(end,k-1))
                continue
            end
% display('existing solutions')
% [kk,k,biomass(end,k)]
% upper and lowber bounds of sum_i t_i
            t_ineq =[repmat(sparse_con(k),NUM,1);repmat(sparse_EX(kk),NUM,1)];
            if NUM==1
                A_ineq = [A_ineq_1; T_ineq];
                b_ineq = [b1_ineq;t_ineq];
            else
                A_ineq = [A_ineq_1; T_ineq; T_ineq2];
                b_ineq = [b1_ineq;  t_ineq; t_ineq2];
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
            flux_all(n+1:n+N_u,k)=1; %LB(n+1:n+N_u)=1;
            model.start =flux_all(:,k);
            %% exact solution
            if EXACT==1
                params.outputflag = 1;%0;
                % params.DisplayInterval=100;
                if (NUM==1 || NUM==3)
                    params.IntFeasTol=1e-9;%Default value:	1e-5
                    params.FeasibilityTol=1e-9;%Default value:	1e-6
    %                 params.MIPGap=1e-9;%0;%
                end
                params.threads=N_cpu;
                params.TimeLimit=10;%3;% the longer the closer to the optimal
                result = gurobi(model,params);
                if  isfield(result,'x')
                flux_all(:,k)=result.x;
                end
                % save models to debug
                % gurobi_write(model, 'iJR904_3.lp');
                % gurobi_write(model, 'iJR904_2.lp');
                % gurobi_write(model, 'iJR904_1.lp');
            elseif EXACT==-1 && NUM==3
                %% heuristic solution using NUM-1 full results
                filename22 = ['iJR904_output_2/iJR904_primal2_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];
                flux_all2=load(filename22,'flux_all');
                flux_all2=flux_all2.flux_all;
                clear params;
                % params.outputflag =0;
                %params.outputflag =1
                params.DisplayInterval=100;
                params.IntFeasTol=1e-9;%Default value:	1e-5
                params.FeasibilityTol=1e-9;%Default value:	1e-6
                % params.MIPGap=0;%1e-9;
                params.TimeLimit=100;
                result =MIP_game0(NUM,UB,model,params,flux_all,flux_all2,k,n,N_u,N_i,N_e);
                if  isfield(result,'x')
                    flux_all(:,k)=result.x;
                end
            else
              %% Similarity-based heuristic solutions/random projection
                for  times=1:21 %repeat%times=1                    
                    mod7=mod(times,7);
                    mod3=mod(times,3);
                    if  (mod3==2) && (mod7==3)
                        continue
                    end
                    if  (mod3==3) && (kk==38)
                        continue
                    end
                    if  (mod3==1) && (kk==1)
                        continue
                    end
                    if mod3==1
                        ex_in_flux=abs(flux_all(N_u+1:n,k))+abs(pre.flux_all(N_u+1:n,min(length(sparse_con),max(k-2+mod7,1))));%
                    elseif mod3==2
                        ex_in_flux=abs(flux_all(N_u+1:n,k))+abs(flux_all(N_u+1:n,min(length(sparse_con),max(k-3+mod7,1))));%
                    else %mod3==3
                        ex_in_flux=abs(flux_all(N_u+1:n,k))+abs(next.flux_all(N_u+1:n,max(k-mod7,1)));%
                    end
                    if biomass(end,k)>0
                        model.start =flux_all(:,k);
                        howmany=150;%the larger the closer to the optimal but slower
                    else
                        model.start =flux_all(:,k-1);
                        howmany=600;%100;
                    end
                    %% inactive intracellular and transport reactions heuristic
                    inactive_ex_in=find(sum(abs(ex_in_flux),2)<10^-7);
                    UB3= UB;
                    UB3(n+N_u+inactive_ex_in)=0;UB3(n+1:n+N_u)=1;
                    UB3(n+N_u+ datasample(inactive_ex_in,howmany,'Replace',false))=1;
                    model.ub = UB3;
                    clear params;
                    params.outputflag = 0;%1%
                    params.DisplayInterval=100;
                    params.IntFeasTol=1e-9;%Default value:	1e-5
                    params.FeasibilityTol=1e-9;%Default value:	1e-6
                    params.MIPGap=0;%1e-9;
                    if biomass(end,k)>0
                        params.TimeLimit=30;
                    else
                    params.TimeLimit=30; % 10;%the larger the closer to the optimal but slower
                    end
                    params.threads=N_cpu;
                    result = gurobi(model,params);
%                     %% Game theory-based heuristic solutions after Similarity-based heuristic
                    if  isfield(result,'x') && mod(times,10)==0%NUM>1 
                        params.TimeLimit=10;
                        result =MIP_game(NUM,UB,result,biomass,model,params,flux_all,k,n,N_u,flux_all_row_id);
                    end
                    %% double check in case of numerical issues
                    if  isfield(result,'x') && biomass(end,k)<result.objval
                        t=result.x(n+1:end);
                        t(1:N_u)=1;
                        cvx_clear
                        %cvx_solver Gurobi_2
                        cvx_precision best %default%low%high%medium%
                        cvx_begin  quiet
                        variable x(n)
                        minimize( norm(x(N_u+1:end),1) )
                        subject to
                        Su*x == b_e(1:end-NUM+1);
                        x_lb2 <= x <= x_ub2;
                        c'*x == sum(result.x(obj>0));
                        diag(x_lb2)*t<= x <= diag(x_ub2)*t;
                        for kkk=1:NUM
                            x(biomass_id(1))==x(biomass_id(kkk));%equal biomass
                        end
                        cvx_end
                        %%
                        if  strcmp( cvx_status,'Solved')%~(isinf(cvx_optval) || isnan(cvx_optval) || abs(cvx_optval)>1000)
                            x =result.x;
                            flux_all(:,k)=result.x;
                            if NUM==1
                                biomass=flux_all(obj>0,:);%
                            else
                                biomass=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
                            end
                            disp( [ kk,k,biomass(end,k)])%,biomass_all_old(kk,k)
                            
                            flux_norm(:,k)=[norm(x(1+N_u:end),1);norm(x,1);norm(x(1:N_u),1)];
                            if save_output1==1
                                save(filename,'sparse_EX','sparse_con','flux_all','biomass','k','kk')
                                disp('saved!!!')
                            end
                        end
                    else
                        continue;
                    end
                end
            end
        end
        if NUM==1
            biomass=flux_all(obj>0,:);
        else
            biomass=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
        end
        clear model;
        for k=1:NUM
            model{k}.sparse_con = sparse_con;
            model{k}.rxns = rxns;
            model{k}.biomass = biomass(k,:);
            model{k}.flux = flux_all(ismember(flux_all_row_id,{'U',['E' num2str( k)],['I' num2str( k)]}),:);
            model{k}.int =  flux_all(ismember(flux_all_row_id,{'TU',['TE' num2str( k)],['TI' num2str( k)]}),:);
        end
        if save_output1==1
        save(filename,'model','sparse_EX','sparse_con','flux_all','rxns','mets','biomass','biomass_id','flux_all_row_id')
        end
        if biomass(1,1)==0
            break
        end
    end
    if save_output2==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use monetone property to update the solutions to speed up
        biomass_all=zeros(length(sparse_EX),length(sparse_con) );
        for kk=1:length(sparse_EX)-1%1
            filename = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];
            load(filename,'biomass')
            biomass_all(kk,1:length(biomass(end,:)))=biomass(end,:);
        end
        %     save('biomass_all3','biomass_all','sparse_con','sparse_EX','NUM')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:length(sparse_EX)-2%1%range_sparseEX%
            filename = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(i) ) '.mat'];
            load(filename)
            for j=1:length(sparse_con)-1
                A=biomass_all(i:end,j:end);
                [M,I] = max(A(:));
                [I_row, I_col] = ind2sub(size(A),I);
                if biomass_all(i,j)<M
                    filename2 = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(i+I_row-1) ) '.mat'];
                    next=load(filename2);
                    flux_all(:,j)=next.flux_all(:,j+I_col-1);
                end
            end
            biomass=[flux_all(biomass_id,:);sum(flux_all(biomass_id,:))];
            clear model;
            for k=1:NUM
                model{k}.sparse_con = sparse_con;
                model{k}.rxns = rxns;
                model{k}.biomass = biomass(k,:);
                model{k}.flux = flux_all(ismember(flux_all_row_id,{'U',['E' num2str( k)],['I' num2str( k)]}),:);
                model{k}.int =  flux_all(ismember(flux_all_row_id,{'TU',['TE' num2str( k)],['TI' num2str( k)]}),:);
            end
            save(filename,'model','sparse_EX','sparse_con','flux_all','rxns','mets','biomass','biomass_id','flux_all_row_id')%'flux_norm',
        end
        toc
    end
end

end