% load original measrured data
fw.model_phantom=input('Please choose a kind of breast phantom (1~999):');
load(['..\data\model' num2str(fw.model_phantom) '\original\Measured_data_0.mat']);
mag0 = Mag_measured_data;
pha0 = Pha_measured_data;
clear Mag_measured_data Pha_measured_data frequency
load(['..\data\model' num2str(fw.model_phantom) '\original\Measured_data_1.mat']);
mag1 = Mag_measured_data;
pha1 = Pha_measured_data;
clear Mag_measured_data Pha_measured_data frequency
background_path = ['..\data\model' num2str(fw.model_phantom) '\background.mat'];
% Calibration_fun(mag0,pha0,mag1,pha1,background_path,freq);
load(['..\data\model' num2str(fw.model_phantom) '\all_material.mat']);
freq = frequency;
Calibration_fun_normalization(mag0, pha0, mag1, pha1, background_path, freq,fw.model_phantom)
