% use GITecoli_iJR904_equal123
clc;
clear;
%%
% number of bacteria
N_ecoli=3%2;%1;%
% save output or not, not needed since output has already been saved
% save_output1 for saving 2nd LP after 1st MIP
% save_output2 for strategy that uses monetone property to update the solutions
save_output1=0;%1;
save_output2=0;%1;
% number of CPU
N_cpu=30%2;
%% EXACT=-1 for NUM==3 GAME heuristic solution from NUM-1==2 results (offer initial solutions for NUM>1)
% EXACT=0 for heuristic but fast, 
% EXACT=1 for exact solution but slow to verify optimality,
EXACT=0%%-1%1;
%%
GITecoli_iJR904_equal123(N_ecoli,save_output1,save_output2,N_cpu,EXACT)