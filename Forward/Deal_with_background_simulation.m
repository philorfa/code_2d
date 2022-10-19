% deal with the original_measurement as background
function Deal_with_background_simulation(fw)

num_file=length(fw.save_path_whole);
for i=1:num_file
    [mag_background,pha_background]=Deal_FDTD_forward(fw.freqs(i)*1e9,fw.save_path_whole{i});
    mag_background_set(:,:,i)=20*log10(mag_background); %dB
    pha_background_set(:,:,i)=pha_background;
end
save(['..\data\model' num2str(fw.model_phantom) '\background.mat'],'mag_background_set','pha_background_set');
disp('The new simulated data have been saved!')
end



