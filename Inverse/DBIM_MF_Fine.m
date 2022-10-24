%==================DBIM for multiple frequency and frequency hopping
close all
clear
% clc
InP.test_name=input('Please enter name of this test:','s');%%%%name folder will be created in the result folder
InP.new_res=input('Please enter a new coarse resolution on a mulitple of 0.5(mm):');%%%resolution
InP.linear_method=input('Please choose a linear method from TwIST, CGLS with L-curve, or LSQR (1, 2, or 3, etc.):');%%%%linear solver. check in inverse solver.m
InP.model_phantom=input('Please choose the model nuber');%%%%for example 146
% maximum total number of iterations
total_freqs{1} = [0.9]*1e9; %frequency hopping 
total_freqs{2} = [1.0]*1e9;
total_freqs{3} = [1.1]*1e9;
total_freqs{4} = [1.3]*1e9;
total_freqs{5} = [1.4]*1e9;

%total_freqs{1} = [0.5 0.8]*1e9; %GHz mulitple frequency 0.5 and 0.8
%total_freqs{1} = [1]*1e9;%%%can be mixed with frequency hopping

maxIter=[10 10 10 10 10];%%number of iterations for each total_freqs index.
% check input, frequency, and path
path = Check_Input(InP, total_freqs);

global multiple_sigma
multiple_sigma = 1;
disp(['Adjustment factor is ' num2str(multiple_sigma)]);
material_flag=[1,1.5];
Result = DBIM_Inverse_Fun(InP,total_freqs,maxIter,'Inverse',1,1,material_flag);
save([ path InP.test_name  '--Result.mat'],'Result');
disp('All results has been saved!');
Save_all_figures(path);

