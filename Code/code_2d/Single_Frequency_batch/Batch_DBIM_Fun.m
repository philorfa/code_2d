function Batch_DBIM_Fun(InP, total_freqs, maxIter,material_flag)

% the same as the code the DBIM_MF_Fine
path = Check_Input(InP, total_freqs);%%%check input value
Result = DBIM_Inverse_Fun(InP,total_freqs,maxIter,'Inverse',1,1,material_flag);%%%%main dbim part
save([ path InP.test_name  '--Result.mat'],'Result');%%%save results
disp('All results has been saved!');
Save_all_figures(path);

end
