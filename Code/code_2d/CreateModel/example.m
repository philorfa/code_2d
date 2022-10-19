%%%%please change to the current file folder
modelname='model10000';
if ~exist(modelname,'dir')==1
   mkdir(modelname);
end

size=500;
ant_rad=65;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=zeros(size);
Debye_fine_test_model_Sigma_s=zeros(size);
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%      =ellipse([size/2,size/2],100,100,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
%  [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%      =ellipse([size/2,size/2],148,190,2,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
[Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
     =ellipse([size/2,size/2],158,218,1,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);

[Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
     =ellipse([size/2,size/2],135,180,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
% for i=1:size
%     for j=1:size
%         if (i-size/2)^2+(j-size/2)^2<=200^2 &&(i-size/2)^2+(j-size/2)^2>=197^2
%         Debye_fine_test_model(i,j)=2;
%         end
%     end
% end
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%     =ellipse([size/2,size/2],160,160,2,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%     =ellipse([size/2,size/2],115,160,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%     =ellipse([size/2,size/2],125,170,4,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%      =ellipse([size/2,size/2],115,160,3,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
% [Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s]...
%     =ellipse([size/2,size/2],115-6,160-6,0.5,Debye_fine_test_model,Debye_fine_test_model_EpsInf,Debye_fine_test_model_EpsDelta,Debye_fine_test_model_Sigma_s,0,[0,0,0]);
for i=1:8
    antLocations_fine_2D(i,1)=size/2+2*ant_rad*(cos(-(i-1)*pi/4+pi/2));
    antLocations_fine_2D(i,2)=size/2+2*ant_rad*(sin(-(i-1)*pi/4+pi/2));
end
a=168;
b=212;
size1=480;
r1=sqrt(2*a^2*b^2/(a^2+b^2));
antLocations_fine_2D=zeros(8,2);
for i=1:2:8
    antLocations_fine_2D(i,1)=size1/2+cos(-2*pi/8*(i-1)+pi/2)*a;
    antLocations_fine_2D(i,2)=size1/2+sin(-2*pi/8*(i-1)+pi/2)*b;
end
for i=2:2:8
    antLocations_fine_2D(i,1)=size1/2+r1*cos(-2*pi/8*(i-1)+pi/2);
    antLocations_fine_2D(i,2)=size1/2+r1*sin(-2*pi/8*(i-1)+pi/2);
end


for i=[1,7]
    antLocations_fine_2D(i,1)=size/2+220*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+222*(sin(-(i-1)*pi/6+pi/2));
end
for i=[2,8]
    antLocations_fine_2D(i,1)=size/2+215*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+215*(sin(-(i-1)*pi/6+pi/2));
end
for i=[3,9]
    antLocations_fine_2D(i,1)=size/2+189*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+189*(sin(-(i-1)*pi/6+pi/2));
end
for i=[4,10]
    antLocations_fine_2D(i,1)=size/2+176*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+176*(sin(-(i-1)*pi/6+pi/2));
end
for i=[5,11]
    antLocations_fine_2D(i,1)=size/2+189*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+189*(sin(-(i-1)*pi/6+pi/2));
end
for i=[6,12]
    antLocations_fine_2D(i,1)=size/2+215*(cos(-(i-1)*pi/6+pi/2));
    antLocations_fine_2D(i,2)=size/2+215*(sin(-(i-1)*pi/6+pi/2));
end
% antLocations_fine_2D(1,1)=size/2;
% antLocations_fine_2D(1,2)=size/2+148;
% antLocations_fine_2D(5,1)=size/2;
% antLocations_fine_2D(5,2)=size/2-148;
% antLocations_fine_2D(2,1)=size/2+130/sqrt(2);
% antLocations_fine_2D(2,2)=size/2+130/sqrt(2);
% antLocations_fine_2D(6,1)=size/2-130/sqrt(2);
% antLocations_fine_2D(6,2)=size/2-130/sqrt(2);
% antLocations_fine_2D(3,1)=size/2+110;
% antLocations_fine_2D(3,2)=size/2;
% antLocations_fine_2D(7,1)=size/2-110;
% antLocations_fine_2D(7,2)=size/2;
% antLocations_fine_2D(4,1)=size/2+130/sqrt(2);
% antLocations_fine_2D(4,2)=size/2-130/sqrt(2);
% antLocations_fine_2D(8,1)=size/2-130/sqrt(2);
% antLocations_fine_2D(8,2)=size/2+130/sqrt(2);
% antLocations_fine_2D=ceil(antLocations_fine_2D);
imagesc(Debye_fine_test_model');
axis xy
axis square
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2),'r*')

%%%%%%% construct all_material.mat
frequency=(0.5:0.1:2.5)*1e9;
mats=struct;
mats=AddMaterials(mats,'immer',0,[6.566,16.86,0.3231,1.4288e-10]);
%mats=AddMaterials(mats,'immer1',1,[6.566,16.86,0.3231,1.4288e-10]);
mats=AddMaterials(mats,'gly',1,[6.566,16.86,0.3231,1.4288e-10]);
mats=AddMaterials(mats,'brain',0.5,[20,20,0.147,0]);
mats=AddMaterials(mats,'csf',3,[18.28,44.71,0.07822,0]);
mats=AddMaterials(mats,'plastic',2,[3,0,0,0]);

% init_guess=struct;
% init_guess=AddDebye(init_guess,[10,10,0.147,0]);
% save(['./' modelname '/all_material.mat'],'frequency','mats','init_guess')
% save(['./' modelname '/' modelname '.mat'],'Debye_fine_test_model','antLocations_fine_2D','Debye_fine_test_model_EpsInf','Debye_fine_test_model_Sigma_s','Debye_fine_test_model_EpsDelta')
% 
% 
