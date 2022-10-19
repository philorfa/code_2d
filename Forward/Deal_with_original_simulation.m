% deal with the original_measurement as background
function Deal_with_original_simulation(fw)
switch fw.mode
    case 'Forward_original'
        for i=1:length(fw.freqs)
            [mag_background,pha_background]=Deal_FDTD_forward(fw.freqs(i),fw.save_path_whole{1});
            Calibrated_Mag(:,:,i)=20*log10(mag_background); %dB
            Calibrated_Pha(:,:,i)=pha_background; %rad
        end
        save(['..\data\model' num2str(fw.model_phantom) '\Calibrated_data.mat'],'Calibrated_Mag','Calibrated_Pha');
        disp('The new Calibrated data based on the original model have been saved!')
    case 'Forward_background'
        for i=1:length(fw.freqs)
            [mag_background,pha_background]=Deal_FDTD_forward(fw.freqs(i),fw.save_path_whole{1});
            mag_background_set(:,:,i)=20*log10(mag_background); %dB
            pha_background_set(:,:,i)=pha_background; %rad
        end
        save(['..\data\model' num2str(fw.model_phantom) '\background.mat'],'mag_background_set','pha_background_set');
        disp('The new background data have been saved!')
end
end



