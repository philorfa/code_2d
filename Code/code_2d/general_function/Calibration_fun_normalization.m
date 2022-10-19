function Calibration_fun_normalization(mag0, pha0, mag1, pha1, background_path, frequency)
load(background_path);
Diff_Mag=mag1-mag0;
Diff_Pha=pha1-pha0;
[a,b,c]=size(Diff_Pha);
% ratio = 1;

 %ratio_t = min(squeeze(mag_background_set(1,2:end,1:c)))./ min(squeeze(mag0(1,2:end,1:c)));
 %ratio = reshape(repmat(ratio_t,a*b,1),a,b,c);
%Calibrated_Mag=mag1;
%Calibrated_Pha=Pha1;
Calibrated_Mag=Diff_Mag+mag_background_set(:,:,1:c);
Calibrated_Mag(Diff_Mag==0)=0;
%Calibrated_Pha=pha1;
Calibrated_Pha=Diff_Pha+pha_background_set(:,:,1:c);
Calibrated_Pha(Diff_Mag==0)=0;
% output the mag of diff
h=figure;
a=1:size(frequency);
plot(squeeze(Diff_Mag(1,:,1:end)));
grid on
xlabel('Received Ant');
ylabel('Diff (dB)');
Add_legend(h,'freq = ' ,frequency(1:end), ' GHz');
% Add_Marker(h);
save('Calibrated_data','Calibrated_Mag','Calibrated_Pha','frequency');
disp('Calibrated data have been saved in current folder!');
end