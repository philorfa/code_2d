% db_MF_Forward_FDTD
function [save_path_whole] = db_MF_Forward_FDTD(fctr,fw)
tic
load sim_params_2D
%-------------------Forward solution parameter------------------
warning('This is forward solution process !!');
%===========transform the original model into a new resolution============
[antLocations_coarse_2D, resolution, mul, coarse_grid_x, coarse_grid_y,...
    Debye_coarse_2D_model, numAnts] = Transform_Resolution(fw);
%=======temp ant=================
warning(['The number of antenna is ' num2str(numAnts)]);
%=====================initialize parameters=============================
% timeSteps=timeSteps_fine;
timeSteps=6000;
omegas = 2.*pi.*fctr;
%obs_iter = 5;  % used to plot inverse solution
% load dontUse;   % indicates TR pairs NOT to be used in reconstruction
numTR = numAnts*(numAnts-1)/2;
PairTR = Generate_Ant_index(numAnts);
% computational domain dimensions
dimX = coarse_grid_x;%from temp_coarse_2D.mat
dimY = coarse_grid_y;%from temp_coarse_2D.mat
delX = resolution;%from temp_coarse_2D.mat
delT = delX/(2*c); 
%  tauP=1.5e-11;
%========================loading model ==================================
[estEpsInf, estEpsDelta, estCond, tauP, EpsS, bbSize, ~,...
    bbox_interior_mask, ~] = Load_FDTD_material_2D(Debye_coarse_2D_model,fw.model_phantom,fw.mode,0,fw.material_flag);
%===============pre-process the FDTD parameters===================
FDTD_Prep_path = FDTD_Preprocess(num2str(fw.model_phantom), delT,pml_coarse,dimX,dimY,bbox_interior_mask,mu0,eps0,...
    delX,EpsS,antLocations_coarse_2D);
%%%%%%%%%%%%%%%%%DBIM beginning%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start to reconstruct the MWI...')
numFreqs=1;
% compute frequency-domain quantities for current source
switch fw.mode
    case 'Forward_original'
        [source, timeSteps]=create_GMS(fctr,delT,'wide',timeSteps,dimX,dimY,0);
%           source = create_GMS2( fctr, timeSteps, delT, delay, tauS, mul);
    case 'Forward_background'
        [source, timeSteps]=create_GMS(fctr,delT,'wide',timeSteps,dimX,dimY,0);
end
[esource_mag,esource_pha]=fft_trans(source,omegas,delT,0,0);
disp(['The current frequency is ' num2str(fctr/1e9) 'GHz'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DBIM iterations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============== run function for FDTD simulations ==============
disp( ' > forward solver ...(FDTD is running)')
%save('a.mat','FDTD_Prep_path','estEpsInf', 'estEpsDelta','estCond', 'tauP','delT', 'eps0', 'mu0', 'delX', 'source', 'numAnts', 'numFreqs','omegas','timeSteps', 'EpsS', 'bbSize', 'esource_mag', 'esource_pha');
[greensFunctions_mag, greensFunctions_pha, ~,~, original_measurement] = project_FDTD_2D(FDTD_Prep_path, estEpsInf, estEpsDelta, ...
    estCond, tauP, delT, eps0, mu0, delX, source, numAnts, numFreqs, ...
    omegas,timeSteps, EpsS, bbSize, esource_mag, esource_pha);

%save(num2str(fw.fctr/1e8),'greensFunctions_mag', 'greensFunctions_pha');
%--------------return for the foward process----------------
disp('The measured data has been saved.');
save_path=['..\data\model' num2str(fw.model_phantom) '\model'...
    num2str(fw.model_phantom) fw.mode];
mkdir(save_path);
save_path_whole=[save_path '\simulated_model' num2str(fw.model_phantom) '_' ...
    num2str(fctr/1e9) 'GHz.mat'];
save(save_path_whole,'original_measurement', 'delT','delX','source');
delete(FDTD_Prep_path);
toc
end
