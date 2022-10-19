clear
size=200;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=ones(size);
Debye_fine_test_model_Sigma_s=zeros(size);
for i=60:140
    for j=60:140
        %if (i-size/2)^2+(j-size/2)^2<=40^2
            Debye_fine_test_model(i,j)=0.5;
        %end
    end
end
            
for i=80:100
    for j=100:120
        Debye_fine_test_model(i,j)=0.5;
        Debye_fine_test_model_EpsInf(i,j)=1.6;
        %Debye_fine_test_model_Sigma_s(i,j)=0.1;

    end
end
for i=110:130
    for j=100:120
        Debye_fine_test_model(i,j)=0.5;
        Debye_fine_test_model_EpsInf(i,j)=1.6;
        %Debye_fine_test_model_Sigma_s(i,j)=0.1;
    end
end
for i=1:8
    antLocations_fine_2D(i,1)=400+4*75*(cos((i-1)*pi/4));
    antLocations_fine_2D(i,2)=400+4*75*(sin((i-1)*pi/4));
end
Debye_fine_test_model=imresize(Debye_fine_test_model,4,'nearest');
Debye_fine_test_model_EpsInf=imresize(Debye_fine_test_model_EpsInf,4,'nearest');
Debye_fine_test_model_EpsDelta=imresize(Debye_fine_test_model_EpsDelta,4,'nearest');
Debye_fine_test_model_Sigma_s=imresize(Debye_fine_test_model_Sigma_s,4,'nearest');
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