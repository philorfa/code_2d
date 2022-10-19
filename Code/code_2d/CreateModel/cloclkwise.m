%%%clockwise circle array model

size=480;
Debye_fine_test_model=zeros(size);
Debye_fine_test_model_EpsDelta=zeros(size);
Debye_fine_test_model_EpsInf=zeros(size);
Debye_fine_test_model_Sigma_s=zeros(size);
for i=1:size
    for j=1:size
        if (i-size/2)^2/170^2+(j-size/2)^2/170^2<=1
            Debye_fine_test_model(i,j)=0.5;
        end
    end
end
for i=1:8
    antLocations_fine_2D(i,1)=ceil(size/2+cos(-2*pi/8*(i-1)+pi/2)*180);
    antLocations_fine_2D(i,2)=ceil(size/2+sin(-2*pi/8*(i-1)+pi/2)*180);
end
imagesc(Debye_fine_test_model');
axis xy
axis square
hold on
scatter(antLocations_fine_2D(:,1),antLocations_fine_2D(:,2))