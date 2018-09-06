% use GITecoli_core_equal12
clc;
clear;
%%
% number of bacteria
N_ecoli=2;%1;%
% save output or not, not needed since output has already been saved
% save_output1 for 1st MIP
% save_output2 for 2nd LP
save_output1=0;%1;
save_output2=0;%1;
% number of CPU
N_cpu=2;
%%
GITecoli_core_equal12(N_ecoli,save_output1,save_output2, N_cpu)