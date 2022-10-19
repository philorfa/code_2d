%%please change to the current file folder
modelname='model11111';
if ~exist(modelname,'dir')==1
   mkdir(modelname);
end

size=400;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=zeros(size);
Debye_fine_test_model_Sigma_s=zeros(size);


[Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
    =ellipse([size/2,size/2],100,100,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);

for i=1:8
    antLocations_fine_2D(i,1)=size/2+66*2*cos((i-1)*pi/4+pi/2);
    antLocations_fine_2D(i,2)=size/2+66*2*sin((i-1)*pi/4+pi/2);
end
    imagesc(Debye_fine_test_model');
axis xy
axis square
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2))

%%%%%%% construct all_material.mat
save(['./' modelname '/all_material.mat'],'frequency','mats','init_guess')
save(['./' modelname '/' modelname '.mat'],'Debye_fine_test_model','antLocations_fine_2D','Debye_fine_test_model_EpsInf','Debye_fine_test_model_Sigma_s','Debye_fine_test_model_EpsDelta')


