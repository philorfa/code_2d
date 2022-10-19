%%%%please change to the current file folder
modelname='model10000';
if ~exist(modelname,'dir')==1
   mkdir(modelname);
end

size=480;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=zeros(size);
Debye_fine_test_model_Sigma_s=zeros(size);


[Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
    =ellipse([size/2,size/2],110,148,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);

for i=1:8
antLocations_fine_2D(i,1)=size/2+180*cos(-(i-1)*pi/4+pi/2);
antLocations_fine_2D(i,2)=size/2+180*sin(-(i-1)*pi/4+pi/2);
end
antLocations_fine_2D=ceil(antLocations_fine_2D);
imagesc(Debye_fine_test_model');
axis xy
axis square
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2))

%%%%%%% construct all_material.mat
frequency=(0.5:0.1:2.5)*1e9;
mats=struct;
mats=AddMaterials(mats,'immer',0,[6.566,16.86,0.3231,1.4288e-10]);
mats=AddMaterials(mats,'brain',0.5,[21,20,0.147,0]);

init_guess=struct;
init_guess=AddDebye(init_guess,[21,20,0.147,0]);
save(['./' modelname '/all_material.mat'],'frequency','mats','init_guess')
save(['./' modelname '/' modelname '.mat'],'Debye_fine_test_model','antLocations_fine_2D','Debye_fine_test_model_EpsInf','Debye_fine_test_model_Sigma_s','Debye_fine_test_model_EpsDelta')


