%DBIM_MF_Fine
close all
clear
clc
dbstop if all error
InP.test_name=input('Please enter name of this test:','s');
InP.new_res=input('Please enter a new coarse resolution on a mulitple of 0.5(mm):');
InP.logic_skin=input('Please choose if skin will be known or not( 1 or 0 ):');
InP.opt_init=input('Please choose whether to optimize initial guess(1 or 0):');
% InP.SNR=input('Please enter the SNR of signal(0~150,1000 for no noise):');
InP.Opt_L1=input('Please decide whether to activate the L1-norm regularization (1 or 0):');
InP.linear_method=input('Please choose a linear method from TwIST and CGLS (1 or 2):');
InP.model_phantom=input('Please choose a kind of breast phantom ( 1~99 ):');
% maximum total number of iterations
total_freqs{1} = [2.0]*1e9;
%total_freqs{2} = [1.6]*1e9;
%total_freqs{3} = [1.9]*1e9;
%total_freqs{4} = [2.3]*1e9;
%total_freqs{5} = [2.6]*1e9;
%total_freqs{6} = [2.9]*1e9;
 %total_freqs{4} = [21]*1e9;
% total_freqs{5} = [26]*1e9;
% total_freqs{6}=[16]*1e9;
maxIter=[15];
% check input, frequency, and path
path = Check_Input(InP, total_freqs);


switch InP.opt_init
    case 0 % start DBIM without optimized initial guess
        Result = DBIM_Inverse_Funnodispersive(InP,total_freqs,maxIter,'Inverse',1,1);
        save([ path InP.test_name  '--Result.mat'],'Result');
        disp('All results has been saved!');
        Save_all_figures(path);
    case 1 %start DBIM with optimized initial guess
        maxIter_opt=input('Please input the number of iterations for changing the initial guess:');
        load(['..\data\model' num2str(InP.model_phantom) '\all_material.mat']);
        num_sample=length(init_guess.EpsInf);
        for i=1:num_sample
            init_result(i) = DBIM_Inverse_Funnodispersive(InP,total_freqs(1),maxIter_opt,'Inverse',i,0);
        end
        plot_initial_guess(init_result,init_guess);
        opt_init = Analyze_optimized_init(init_result,init_guess);
        save([ path InP.test_name '--init_result.mat'],'init_result','opt_init');
        Save_all_figures(path);
        close all
        Result = DBIM_Inverse_Funnodispersive(InP,total_freqs,maxIter,'Inverse',opt_init,1);
        save([ path InP.test_name '--Result.mat'],'Result');
        disp('All results has been saved!');
        Save_all_figures(path);
end
