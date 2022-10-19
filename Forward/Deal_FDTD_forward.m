function [mag_background,pha_background] = Deal_FDTD_forward(fctr,path_t)
load(path_t);
[numAnts,~,len] = size(original_measurement);
omega = 2*pi*fctr;

[esource_mag,esource_pha] = fft_trans(source,omega,delT,0,0);
[mag,pha] = fft_trans(reshape(original_measurement,numAnts^2,len),...
   omega,delT,esource_mag,esource_pha);
mag_background = reshape(mag,numAnts,numAnts);
pha_background = reshape(pha,numAnts,numAnts);
end
