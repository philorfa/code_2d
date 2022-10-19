function  [antLocations_coarse_2D, resolution, mul, coarse_grid_x, coarse_grid_y,...
    Debye_coarse_2D_model,numAnts] = Transform_Resolution(input_struct)
% input_struct may be InP or fw.
% The function is important to initize the coarse model from fine model.
% All coarse data are saved in the file of temp_coarse_2D.
def_res=0.5; %mm
if (mod(input_struct.new_res,def_res)~=0)
    disp('resolution is not integral multiple of 0.5 mm')
    return
end
mul=input_struct.new_res/def_res;%%%%multiplifier

load(['..\data\model' num2str(input_struct.model_phantom) '\model' num2str(input_struct.model_phantom) '.mat'])%%%%load mmodel 
antLocations_coarse_2D=ceil(antLocations_fine_2D/mul);
[numAnts,~]=size(antLocations_coarse_2D);

Debye_coarse_2D_model.model=imresize(Debye_fine_test_model,1/mul,'nearest');
Debye_coarse_2D_model.DeltaEps=imresize(Debye_fine_test_model_EpsDelta,1/mul,'nearest');
Debye_coarse_2D_model.EpsInf=imresize(Debye_fine_test_model_EpsInf,1/mul,'nearest');
Debye_coarse_2D_model.SigmaS=imresize(Debye_fine_test_model_Sigma_s,1/mul,'nearest');

resolution=input_struct.new_res/10^3;
[coarse_grid_x,coarse_grid_y]=size(Debye_coarse_2D_model.model);

end
