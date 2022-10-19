%%%%please change to the current file folder
modelname='model10000';
if ~exist(modelname,'dir')==1
   mkdir(modelname);
end
%%%%%%% construct all_material.mat
frequency=(0.5:0.1:3.5)*1e9;
mats=struct;
mats=AddMaterials(mats,'immer',0,[2.600000000000000,0,0,1.712500000000000e-11]);
mats=AddMaterials(mats,'skin',2,[17.9872,22.2744,0.738300000000000,1.712500000000000e-11]);
init_guess=struct;
init_guess=AddDebye(init_guess,[5.76,5.51,0.080200000000000,1.712500000000000e-11]);

save(['./' modelname '/all_material.mat'],'frequency','mats','init_guess')
save(['./' modelname '/' modelname '.mat'],'Debye_fine_test_model','antLocations_fine_2D','Debye_fine_test_model_EpsInf','Debye_fine_test_model_Sigma_s','Debye_fine_test_model_EpsDelta')


