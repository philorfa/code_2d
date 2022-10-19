% load original measrured data
load ('./original/Measured_data_0.mat');
mag0 = Mag_measured_data;
pha0 = Pha_measured_data;

freq = frequency;
clear Mag_measured_data Pha_measured_data frequency
load ('./original/Measured_data_1.mat');
mag1 = Mag_measured_data;
pha1 = Pha_measured_data;
clear Mag_measured_data Pha_measured_data frequency
background_path = './background.mat';
% Calibration_fun(mag0,pha0,mag1,pha1,background_path,freq);
Calibration_fun_normalization(mag0, pha0, mag1, pha1, background_path, freq)
