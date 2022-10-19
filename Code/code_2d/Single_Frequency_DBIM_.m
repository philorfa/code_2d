% function Single_Frequency_DBIM_2D
clear
close all
time_total=tic;%%%time rec start
freq_test=[0.5:0.1:3.0]; %%%frequency range (same as the one in the model folder)
maxIter=20;%%%number of iteration for each frequency 
 
for ind=freq_test(6)%%%% frequency index for reconstruction

%=================input================
total_freqs{1}=ind*1e9;
InP.model_phantom=26;%%%model index
k=find(freq_test==ind);
InP.test_name=['test' num2str(InP.model_phantom*100+k)];
InP.new_res=2;%%%resolution
InP.linear_method=6;%%%%select inverse solver
material_flag=[1];%%%%define flags of the reconstruction area (flag are defined in the model***.mat)

Batch_DBIM_Fun(InP, total_freqs, maxIter,material_flag);%%main DBIM 
close all
end
time_t=toc(time_total);%%%%time record end
path=['..\result\' InP.test_name '\'];%%%save path for the time
save([ path InP.test_name '--time.mat'],'time_t');