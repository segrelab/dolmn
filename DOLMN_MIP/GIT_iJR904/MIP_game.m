function result =MIP_game(NUM,UB,result,biomass,model,params,flux_all,k,n,N_u,flux_all_row_id)
%% game theory heuristic MIP_game.m
% solve sub-MIP which is small and simple
%%
UB3= UB;
if  isfield(result,'x') && biomass(end,k)<result.objval
    model.start =result.x;
    temp=result.x(n+N_u+1:n+n);
else
    temp=flux_all(n+N_u+1:n+n,k);
end
UB3(n+N_u+1:n+n)=temp;
UB3(ismember(flux_all_row_id,{'TI1'}))=1;
model.ub = UB3;
result = gurobi(model,params);
%%
UB3= UB;
if  isfield(result,'x') && biomass(end,k)<result.objval
    model.start =result.x;
    temp=result.x(n+N_u+1:n+n);
else
    temp=flux_all(n+N_u+1:n+n,k);
end
UB3(n+N_u+1:n+n)=temp;
UB3(ismember(flux_all_row_id,{'TE1'}))=1;
model.ub = UB3;
result = gurobi(model,params);

if NUM>=2
    %%
    UB3= UB;
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TE1','TE2'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
    %%
    UB3= UB;
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TI1','TE1'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
    %%
    UB3= UB;
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TI2','TE2'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
end
if NUM>=3
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
    %%
    UB3= UB;
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TE1','TE3'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
    %%
    UB3= UB;
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TE2','TE3'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
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
    %%
    UB3= UB;%sum(UB(n+1:n+N_u))
    if  isfield(result,'x') && biomass(end,k)<result.objval
        model.start =result.x;
        temp=result.x(n+N_u+1:n+n);
    else
        temp=flux_all(n+N_u+1:n+n,k);
    end
    UB3(n+N_u+1:n+n)=temp;
    UB3(ismember(flux_all_row_id,{'TI3','TE3'}))=1;
    model.ub = UB3;
    result = gurobi(model,params);
end
end
