%%%FDTD forward code
clear
fw.new_res=input('Please enter a new coarse resolution on a mulitple of 0.5(mm):');
fw.model_phantom=input('Please choose a kind of breast phantom (1~999):');
fw.mode=input('Please input the mode of forward FDTD: \n 1.Known_original \n 2.Unknown_original \n');
load(['..\data\model' num2str(fw.model_phantom) '\all_material.mat'],'frequency');
fw.freqs=frequency;
switch fw.mode
    case 1 %%%%when FDTD data is used as measurement
        fw.fctr=1.0e9; % 2.0 GHZ
        fw.mode='Forward_original';
        fw.material_flag=[0];
        disp(['The central frequency is ' num2str(fw.fctr/1e9) ' GHz']);
        fw.save_path_whole{1} = db_MF_Forward_FDTD(fw.fctr,fw);
        disp('Starting to deal with the simulated data...');
        Deal_with_original_simulation(fw); 
    case 2%%%%used  to calculate background.mat, i.e. forward solver.
        % wide band model
        fw.fctr=1e9; % 2.0 GHZ
        fw.mode='Forward_background';
        fw.material_flag=[0,1,2,3,4,5];
        disp(['The central frequency is ' num2str(fw.fctr/1e9) ' GHz']);
        fw.save_path_whole{1} = db_MF_Forward_FDTD(fw.fctr,fw);
        disp('Starting to deal with the simulated data...');
        Deal_with_original_simulation(fw);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % narrow band model
% % %         fw.mode='Forward_background'; %without original model
% % %         for ind=1:length(fw.freqs)
% % %             fw.save_path_whole{ind} = db_MF_Forward_FDTD(fw.freqs(ind),fw);
% % %             disp(fw.save_path_whole{ind});
% % %         end
% % %         disp('Starting to deal with the simulated data...');
% % %         Deal_with_background_simulation(fw); % used for befor we try to reconstruct from External measured data
    otherwise
        error('Please check the input of the forward mode');
end

save(['..\data\model' num2str(fw.model_phantom) '\forward_para.mat'],'fw');
disp('forward parameters have been saved!');