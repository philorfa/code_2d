clear
size=2400;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=ones(size);
Debye_fine_test_model_Sigma_s=zeros(size);
for i=600:1800
    for j=600:1800
        %if (i-size/2)^2+(j-size/2)^2<=40^2
            Debye_fine_test_model(i,j)=1.5;
        %end
    end
end
            
for i=1000:1400
    for j=1000:1400
        if (i-1200)^2+(j-1200)^2<=90^2
        Debye_fine_test_model(i,j)=0.5;
        Debye_fine_test_model_EpsInf(i,j)=1.6;
        %Debye_fine_test_model_Sigma_s(i,j)=0.1;
        end
    end
end
for i=1000:1600
    for j=1000:1600
        if (i-1440)^2+(j-1200)^2<=90^2
        Debye_fine_test_model(i,j)=0.5;
        Debye_fine_test_model_EpsInf(i,j)=1.6;
        %Debye_fine_test_model_Sigma_s(i,j)=0.1;
        end
    end
end
for i=1:8
    antLocations_fine_2D(i,1)=1200+900*(cos((i-1)*pi/4));
    antLocations_fine_2D(i,2)=1200+900*(sin((i-1)*pi/4));
end
Debye_fine_test_model=imresize(Debye_fine_test_model,1,'nearest');
Debye_fine_test_model_EpsDelta=imresize(Debye_fine_test_model_EpsDelta,1,'nearest');
Debye_fine_test_model_EpsInf=imresize(Debye_fine_test_model_EpsInf,1,'nearest');
Debye_fine_test_model_Sigma_s=imresize(Debye_fine_test_model_Sigma_s,1,'nearest');
antLocations_fine_2D=ceil(antLocations_fine_2D*1);
imagesc(Debye_fine_test_model');

axis xy
axis square
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2))


frequency=(0.5:0.1:2.5)*1e9;
mats=struct;
mats=AddMaterials(mats,'immer',0,[1,0,0,1.4288e-10]);
init_guess=struct;
init_guess=AddDebye(init_guess,[1,0,0,0]);