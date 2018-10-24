function result =MIP_game0(NUM,UB,model,params,flux_all,flux_all2,k,n,N_u,N_i,N_e)%result,biomass,,flux_all_row_id
%% game theory heuristic MIP_game0.m
% solve sub-MIP which is small and simple
UB(n+1:n+N_u)=1;
UB3= UB;
if k>1   
    UB3(n+1:n+size(flux_all2,1)/2)=abs(flux_all2(1:size(flux_all2,1)/2,k))>10^-7;
    model.ub = UB3;
    model.start =flux_all(:,k);
    result = gurobi(model,params);
else
    %only fix 1st N_i
    UB3(n+N_u+N_e+1:n+N_u+N_e+N_i)=abs(flux_all2(1+N_u+N_e:N_u+N_e+N_i,k))>10^-7;
    model.ub = UB3;
    model.start =flux_all(:,k);
    result = gurobi(model,params);
    if  ~isfield(result,'x')
        %only fix 2nd N_i
        UB3= UB;
        UB3(n+N_u+N_e+N_i+N_e+1:n+N_u+N_e+N_i+N_e+N_i)=abs(flux_all2(1+N_u+N_e+N_i+N_e:N_u+N_e+N_i+N_e+N_i,k))>10^-7;
        model.ub = UB3;
        model.start =flux_all(:,k);
        result = gurobi(model,params);
    end
end


% if NUM>=2
%     %%
%     UB3= UB;
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TE1','TE2'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
%     %%
%     UB3= UB;
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TI1','TE1'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
%     %%
%     UB3= UB;
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TI2','TE2'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
% end
% if NUM>=3
% %     UB3= UB;
% %     if  isfield(result,'x') && biomass(end,k)<result.objval
% %         model.start =result.x;
% %         temp=result.x(n+N_u+1:n+n);
% %     else
% %         temp=flux_all(n+N_u+1:n+n,k);
% %     end
% %     UB3(n+N_u+1:n+n)=temp;
% %     UB3(ismember(flux_all_row_id,{'TE1','TE2'}))=1;
% %     model.ub = UB3;
% %     result = gurobi(model,params);
%     %%
%     UB3= UB;
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TE1','TE3'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
%     %%
%     UB3= UB;
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TE2','TE3'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
% %     %%
% %     UB3= UB;
% %     if  isfield(result,'x') && biomass(end,k)<result.objval
% %         model.start =result.x;
% %         temp=result.x(n+N_u+1:n+n);
% %     else
% %         temp=flux_all(n+N_u+1:n+n,k);
% %     end
% %     UB3(n+N_u+1:n+n)=temp;
% %     UB3(ismember(flux_all_row_id,{'TI1','TE1'}))=1;
% %     model.ub = UB3;
% %     result = gurobi(model,params);
% %     %%
% %     UB3= UB;
% %     if  isfield(result,'x') && biomass(end,k)<result.objval
% %         model.start =result.x;
% %         temp=result.x(n+N_u+1:n+n);
% %     else
% %         temp=flux_all(n+N_u+1:n+n,k);
% %     end
% %     UB3(n+N_u+1:n+n)=temp;
% %     UB3(ismember(flux_all_row_id,{'TI2','TE2'}))=1;
% %     model.ub = UB3;
% %     result = gurobi(model,params);
%     %%
%     UB3= UB;%sum(UB(n+1:n+N_u))
%     if  isfield(result,'x') && biomass(end,k)<result.objval
%         model.start =result.x;
%         temp=result.x(n+N_u+1:n+n);
%     else
%         temp=flux_all(n+N_u+1:n+n,k);
%     end
%     UB3(n+N_u+1:n+n)=temp;
%     UB3(ismember(flux_all_row_id,{'TI3','TE3'}))=1;
%     model.ub = UB3;
%     result = gurobi(model,params);
% end
end
