function ratio = Calibration_for_normalization(mag0, mag_background_set)
[a,b,c]=size(mag0);
ratio_t = min(squeeze(mag_background_set(1,2:end,1:c)))./ min(squeeze(mag0(1,2:end,1:c)));
ratio = reshape(repmat(ratio_t,a*b,1),a,b,c);

end