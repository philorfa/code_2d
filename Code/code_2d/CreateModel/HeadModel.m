%%CST head model

size=600;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=zeros(size);
Debye_fine_test_model_Sigma_s=zeros(size);
for i=1:size
    for j=1:size
        if (j-340)^2/(205^2)+(i-size/2)^2/(148^2)<=1
            Debye_fine_test_model(i,j)=3;
        end
        if (j-360)^2/(180^2)+(i-size/2)^2/(120^2)<=1
            Debye_fine_test_model(i,j)=4;

        end
    end
end
for i=1:8
    antLocations_fine_2D(i,1)=size/2+cos(2*pi/8*(i-1)+pi/2)*256;
    antLocations_fine_2D(i,2)=size/2+sin(2*pi/8*(i-1)+pi/2)*256;
end
antLocations_fine_2D=ceil(antLocations_fine_2D);
imagesc(Debye_fine_test_model');
axis square
axis xy
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2))