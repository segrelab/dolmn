%core_CVX_minfluxaftermaxbiomass_primal1

clc;
clear;
% load('../iJR904_newS.mat','newS','mets','rxns','N_u','N_e','N_i','M_i','M_e','lb','ub','c');
load('iJR904_input/iJR904_newS20171006.mat','newS','mets','rxns','N_u','N_e','N_i','M_i','M_e','lb','ub','c');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust bound to speed up
x_lb=lb;% %
x_lb(c==1)=0.1; %find(c==1)% biomass>0.1 for each bacteria
unique_x_lb=unique(x_lb)'; % -100.0000  -10.0000         0    8.3900
id_ATPM=find(ismember(x_lb,x_lb(x_lb> x_lb(c==1) ) ) )
rxns( id_ATPM )% 'ATPM'
%rxns(ismember(x_lb,unique_x_lb(3)) ) % 48 in total
x_ub=ub;% unique(x_ub)
Max = 100%1e2;% all   1000
x_lb(x_lb<-Max)=-Max;
x_ub(x_ub>Max)=Max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sx=0 ==> S(:,-53)x(-53)=-S(:,53)x(53) due to x(53)==8.39
%name also need to change!!
% from 1 S to community S
NUM=3;% 2%number of bacteria
newS_store=newS;% used for b_e

newS(:,id_ATPM)=[];
x_lb(id_ATPM)=[]; x_ub(id_ATPM)=[];
c(id_ATPM)=[];
rxns(id_ATPM)=[];
N_i=N_i-1;

Su=sparse( M_e+NUM*M_i,N_u+NUM*(N_i+N_e) );%model1.model.S;
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
find_c=biomass_id;
%equality constraints
A_e = sparse(m,2*n);
A_e(1:m,1:n) = Su;
temp = full( -newS_store(:,id_ATPM)*lb(id_ATPM));% one ATMP
b_e =[temp(1:M_e); repmat(temp(M_e+1:end),NUM,1)];%b_e =[temp(1:M_e); repmat(temp(N_u+1:end),NUM,1)];


if NUM==1
    sparse_con=[N_i,296:-1:250];%[N_i,32:-1:20];%[N_i, 500, 250];%%[N_i, 1000, 500, 250, 125];%%[400,300];%fliplr(300:50:450)
    sparse_EX=[N_e,38,34,24,22:-1:8];
elseif NUM==2
    sparse_con=[N_i,285:-1:214];%sparse_con=[N_i,278, 270, 241:-1:214];%[N_i,30:-1:15];
    sparse_EX=[N_e,45:-1:8];
else%if NUM==2
    sparse_con=[N_i,285:-1:194];%sparse_con=[N_i,278, 270, 241:-1:214];%[N_i,30:-1:15];
    sparse_EX=[N_e,45:-1:8];%sparse_EX'
end

RERUN_ID=[];
for  kk=1:length(sparse_EX)-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flux_all=zeros(n+n,size(sparse_con,2));
    filename = ['iJR904_output_' num2str( NUM ) '/iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];
    filenameCVX = ['iJR904_output_' num2str( NUM ) '/CVX_iJR904_primal' num2str( NUM ) '_TrsptCon_' num2str( sparse_EX(kk) ) '.mat'];    
    CVX=load(filenameCVX,'flux_all');
    
    opt=load(filename);
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
    flux_all_row_id=cell(size(flux_all,1),1);%repmat({'filename.mat'},1,5)
    flux_all_row_id(1:N_u)={'U'};
    flux_all_row_id(n+1:n+N_u)={'TU'};
    for k=1:NUM
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+1:N_u+(N_e+N_i)*(k-1)+N_e)={['E' num2str( k)]};
        flux_all_row_id(N_u+(N_e+N_i)*(k-1)+N_e+1:N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['I' num2str( k)]};
        %flux_all_row_id(n+(N_i)*(k-1)+1:n+(N_i)*(k-1)+N_i)={['T' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+1:n+N_u+(N_e+N_i)*(k-1)+N_e)={['TE' num2str( k)]};
        flux_all_row_id(n+N_u+(N_e+N_i)*(k-1)+N_e+1:n+N_u+(N_e+N_i)*(k-1)+N_e+N_i)={['TI' num2str( k)]};
    end
    
    for k=1:length(sparse_con)%
        %         if sum(start.flux_all(obj>0,k))==0
        %             break
        %         end
        % upper and lowber bounds of sum_i t_i
        t_a_max = sparse_con(k);%sparse_con(11)        
        %% CVX
        %min flux norm
        if opt.biomass(1,k)>0
            t=opt.flux_all(n+1:end,k);
            t(1:N_u)=1;
            cvx_clear
            %    cvx_solver Gurobi_2% cvx_solver gurobi%which gurobi
            cvx_precision best %default%low%high%medium%
            cvx_begin  quiet
            variable x(n)
            %variable t(n) binary%t(NUM*N_i)
            minimize( norm(x(N_u+1:end),1) )
            subject to
            Su*x == b_e;%temp=abs(Su*x- b_e);temp(temp>10^-8)
            x_lb2 <= x <= x_ub2;%temp=x_lb2 -x;temp(temp>10^-8)   temp=x-x_ub2;temp(temp>10^-8)
            c'*x == opt.biomass(1,k)*NUM;%biomass(k);%biomassMAX(k)-0.01<= c'*x;%
            diag(x_lb2)*t<= x <= diag(x_ub2)*t;
            % sum(( x - diag(x_ub2)*t)>10^-7);
            %temp=x - diag(x_ub2)*t;[find(temp>10^-7),temp(temp>10^-7)]
            %sum(( diag(x_lb2)*t-x)>10^-6);
            %temp=diag(x_lb2)*t-x;[find(temp>10^-7),temp(temp>10^-7)]
            for kkk=1:NUM
                x(find_c(1))==x(find_c(kkk));%equal biomass
            end
            cvx_end
            if  ~strcmp( cvx_status,'Solved')
                cvx_status
            end
            %% if no solution, then
            if  isinf(cvx_optval) || isnan(cvx_optval) || abs(cvx_optval)>1000
                
                disp('2nd try')
                %         [kk,k]
                t(t~=1)=0;
                cvx_clear
                %         cvx_solver_settings('MSK_IPAR_NUM_THREADS',1)
                %         cvx_solver_settings('params.threads',10)
                cvx_precision high%default%best %low%medium%
                cvx_begin  quiet
                variable x(n)
                minimize( norm(x(N_u+1:end),1) )
                subject to
                Su*x == b_e;%temp=abs(Su*x- b_e);temp(temp>10^-8)
                x_lb2 <= x <= x_ub2;%temp=x_lb2 -x;temp(temp>10^-8)   temp=x-x_ub2;temp(temp>10^-8)
                c'*x == opt.biomass(1,k)*NUM;%biomass(k);%biomassMAX(k)-0.01<= c'*x;%
                diag(x_lb2)*t<= x <= diag(x_ub2)*t;
                for kkk=1:NUM
                    x(find_c(1))==x(find_c(kkk));%equal biomass
                end
                cvx_end
            end
            
            if  ~strcmp( cvx_status,'Solved')
                cvx_status
                [kk,k]
                RERUN_ID=[RERUN_ID;[kk,k]];
            end
            %if  ~isnan(cvx_optval)%isfield(result,'x') && biomass(end,k)<result.objval
            %x =result.x;
            flux_all(:,k)=[x;t];
            if NUM==1
                biomass=flux_all(obj>0,:);%
            else
                biomass=[flux_all(obj>0,:);sum(flux_all(obj>0,:))];
            end
            flux_norm(:,k)=[norm(x(1+N_u:end),1);norm(x,1);norm(x(1:N_u),1)];
            
        else
            break%continue;%
        end
    end
    
    clear model;
    for k=1:NUM
        model{k}.sparse_con = sparse_con;
        model{k}.rxns = rxns;%rxns([exch_flux, trspt_flux1, intl_flux1]);
        % [~,~,bio_idx] = intersect(rxns(S_info.bio_idx),model{1}.rxns);
        model{k}.biomass = biomass(k,:);%M1_flux(bio_idx,:);
        model{k}.flux = flux_all(ismember(flux_all_row_id,{'U',['E' num2str( k)],['I' num2str( k)]}),:);%M1_flux;  sum(ismember(flux_all_row_id,{'U','E1','I1'}))
        %model{k}.int =  flux_all(ismember(flux_all_row_id,{['T' num2str( k)]}),:);%M1_int;
        model{k}.int =  flux_all(ismember(flux_all_row_id,{'TU',['TE' num2str( k)],['TI' num2str( k)]}),:);
    end
    flux_all(n+1:n+N_u,:)=abs(flux_all(1:N_u,:)>10^-8);
    %%
    save(filenameCVX,'model','sparse_EX','sparse_con','flux_all','rxns','mets','biomass','biomass_id','flux_norm','flux_all_row_id')
    %     flux_norm    
    if biomass(1,1)==0
        break
    end
end
% save('RERUN_ID_Inaccurate20171016.mat','RERUN_ID')
RERUN_ID
[sparse_EX(RERUN_ID(:,1))',sparse_con(RERUN_ID(:,2))']



