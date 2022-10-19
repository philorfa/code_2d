function Result = DBIM_Inverse_Fun(InP,total_freqs,maxIter,DBIM_mode,init_index,plot_flag)
%  This is a matlab version of code for 2-D inverse scattering imaging
%  written by Zhenzhuang Miao (Supervisor: Dr. Panagiotis Kosmas)
%  Created on 15/05/2017
global multiple_sigma
% Original parameters
load sim_params_2D
% transformed model and relevant parameters.
[antLocations_coarse_2D, resolution, mul, coarse_grid_x, coarse_grid_y,...
    Debye_coarse_2D_model, numAnts] = Transform_Resolution(InP);
%=======disp number of ants=================
warning(['The number of antenna is ' num2str(numAnts)]);
%=====================initialize parameters=============================
% timeSteps = timeSteps_fine; % default timesteps. will be optimized by the function of Create_GMS
timeSteps = 6000;
obs_iter = 5;  % used to plot inverse solution
numTR = numAnts*(numAnts-1)/2;
PairTR = Generate_Ant_index(numAnts);
% computational domain dimensions
dimX = coarse_grid_x;%from temp_coarse_2D.mat
dimY = coarse_grid_y;%from temp_coarse_2D.mat
delX = resolution;%from temp_coarse_2D.mat
delT = delX/(2*c);
%  tauP=1.5e-11;
%=====================Create breast model=============================
[estEpsInf, estEpsDelta, estCond, tauP, EpsS, bbSize, range_c,bbox_interior_mask, constraints]...
    = Load_FDTD_material_2D(Debye_coarse_2D_model, InP.model_phantom, InP.logic_skin, ...
    DBIM_mode,init_index);
%===============pre-process the FDTD parameters===================
FDTD_Prep_path = FDTD_Preprocess(InP.test_name, delT,pml_coarse,dimX,dimY,bbox_interior_mask,...
    mu0,eps0,delX,EpsS, antLocations_coarse_2D);
%=======================load calibrated data========================
load(['..\data\model' num2str(InP.model_phantom) '\Calibrated_data.mat']);
load(['..\data\model' num2str(InP.model_phantom) '\forward_para.mat']);
Calibrated_Mag=10.^(Calibrated_Mag/20);
% check the forward resolution is equal to the inverse resolution, if not, stop
if fw.new_res~= resolution * 1e3
    error('inverse resolution is not the same as forward resolution!');
end
%==================initialization of recorded value======================
Result.residual3  = zeros(sum(maxIter),1);
Result.res_diff_1 = zeros(sum(maxIter),1);
Result.res_diff_2 = zeros(sum(maxIter),1);
Result.res_diff_3 = zeros(sum(maxIter),1);
Result.L_curve_x  = zeros(sum(maxIter),1);
Result.L_curve_Ax_b=zeros(sum(maxIter),1);
contrast           =zeros(bbSize*3,1);
Result.allContrast =zeros(bbSize*3,sum(maxIter));
iter_total=1;  % temp para for counting
%%%%%%%%%%%%%%%%%DBIM beginning%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start to reconstruct the MWI...')
for ind_freq = 1:length(total_freqs)  % frequency hopping approach
    
    fctrs = total_freqs{ind_freq};
    omegas = 2*pi*fctrs;
    numFreqs=length(fctrs);
    disp(['Current frequency is ' num2str(fctrs/1e9) ' GHz']);
    % compute frequency-domain quantities for current source
    switch fw.mode
        case 'Forward_original'
            [source,timeSteps]=create_GMS(fw.fctr,delT,'wide',timeSteps, dimX, dimY, 0);
            %               source = create_GMS2( 2e9, timeSteps, delT, delay, tauS, mul);
        case 'Forward_background'
            [source,timeSteps]=create_GMS(fw.fctr,delT,'wide',timeSteps, dimX, dimY, 0);
    end
    [esource_mag,esource_pha] = fft_trans(source, omegas, delT, 0 ,0);
    %load the calibrated data from forward solution
    [~,freq_index,~]=intersect(round(fw.freqs/1e9,1),round(fctrs/1e9,1));
    measMag=squeeze(Calibrated_Mag(:,:,freq_index));
    measPha=squeeze(Calibrated_Pha(:,:,freq_index));
    disp(['The index of freq is ' num2str(freq_index(:)')]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% DBIM iterations %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter_inner = 1:maxIter(ind_freq)
        % ============== run function for FDTD simulations ==============
        disp(['iteration #' num2str(iter_inner) ' > forward solver ...(FDTD is running)' ])
        %--------------------project_FDTD_2D-----------------------
        [greensFunctions_mag, greensFunctions_pha, receivedFields_mag,...
            receivedFields_pha, ~] = project_FDTD_2D(FDTD_Prep_path, estEpsInf, estEpsDelta, estCond,...
            tauP, delT, eps0, mu0, delX, source, numAnts, numFreqs, omegas,timeSteps, EpsS,...
            bbSize, esource_mag, esource_pha);
        %=================inverse solution setup ==================================
        disp(['iteration #' num2str(iter_inner) ' > inverse solver ...'])
        [ Mat_A, Data_ant] = Generate_Matrix( receivedFields_mag, ...
            receivedFields_pha, greensFunctions_mag, greensFunctions_pha, measMag, ...
            measPha, numFreqs, numTR, omegas, eps0, PairTR, tauP, bbSize);
        %==========save residual and relative error before reconstruction=========
        % norm of residual
        Result.residual3(iter_total) = norm(Data_ant,2);
        % norm of difference from coarse model
        Result.res_diff_1(iter_total) = norm(estEpsInf(bbox_interior_mask)...
            -Debye_coarse_2D_model.EpsInf(bbox_interior_mask))./...
            norm(Debye_coarse_2D_model.EpsInf(bbox_interior_mask));
        Result.res_diff_2(iter_total) = norm(estEpsDelta(bbox_interior_mask)...
            -Debye_coarse_2D_model.DeltaEps(bbox_interior_mask))./...
            norm(Debye_coarse_2D_model.DeltaEps(bbox_interior_mask));
        Result.res_diff_3(iter_total) = norm(estCond(bbox_interior_mask)...
            -Debye_coarse_2D_model.SigmaS(bbox_interior_mask))./...
            norm(Debye_coarse_2D_model.SigmaS(bbox_interior_mask));
        %==============Start the solve the linear equation=================
        Solve_Linear_Equation
        % used for adjustment factor
        if numFreqs~=1
            contrast(2*end/3+1:end)=multiple_sigma*contrast(2*end/3+1:end);
        end
        %--------monitor the structure of the solution of every TWIST--------
        %Result.allContrast(:,iter_total) = contrast;
        contrast(2*end/3+1:end) = contrast(2*end/3+1:end).*2*pi*1e9*eps0;
        
        if isnan(contrast)
            error('The value of contrast includes NaN !');
        end
        
        disp(['Residual error = ' num2str(Result.residual3(iter_inner,1))])
        % ============== Update dielectric properties estimate================
        estEpsInf(bbox_interior_mask) = estEpsInf(bbox_interior_mask) + ...
            contrast(1:end/3);
        estEpsDelta(bbox_interior_mask) = estEpsDelta(bbox_interior_mask) + ...
            contrast(end/3+1:2*end/3);
        estCond(bbox_interior_mask) = estCond(bbox_interior_mask) + ...
            contrast(2*end/3+1:end);
        %%==========================constraints==============================
        %  Make sure the Debye parameters are in reasonable range
        %%%%%%%% this code is only for the real breast model%%%%%%%%%%%
%         switch fw.mode
%             % The first case is only used for 4 real breast phantoms
%             case 'Forward_original'
% %                 warning('Hard constraint is created on real breast phantom')
% %                 estEpsInf((estEpsInf<2.28)&range_c) = 2.28;
% %                 estEpsDelta(estEpsDelta<1.3&range_c) = 1.3;
% %                 estCond(estCond< 0.0023&range_c) = 0.0023;
% %                 estEpsInf(estEpsInf>23.2&range_c) = 23.2;
% %                 estEpsDelta(estEpsDelta>38&range_c) = 38;
% %                 estCond(estCond>0.8&range_c) = 0.8;
%                 estEpsInf(estEpsInf< 3.512 & range_c) = 3.512;
%                 estEpsDelta(estEpsDelta< 0 & range_c) = 0; % for water
%                 estCond(estCond< 0.0655 & range_c) = 0.0655;
%                 estEpsInf(estEpsInf>90&range_c) = 90; % for water
%                 estCond(estCond>1.6&range_c) = 1.6; % for water
%             case 'Forward_background'
%                 % range_c has been created by the function 'Create_model'.
%                 %load(['..\data\model' num2str(InP.model_phantom) '\all_material.mat']);
%                 estEpsInf((estEpsInf < constraints.EpsInf) & range_c) = constraints.EpsInf;
%                 estEpsDelta((estEpsDelta < constraints.DeltaEps) & range_c) = constraints.DeltaEps;
%                 estCond(estCond < constraints.SigmaS & range_c) = constraints.SigmaS;
%                 estEpsInf(estEpsInf > 90 & range_c) = 90;
%         
%%%for cst brain model
% Einf_imm=6.566;
% Edelta_imm=16.86;
% Cond_imm=0.3231;
%                   estEpsInf(estEpsInf< 1& range_c) = 1;
% %                 %estEpsDelta(estEpsDelta> Edelta_imm & range_c) = Edelta_imm; % for material
%                  estCond(estCond~= 0 & range_c) =0;
%                  %estEpsInf(estEpsInf>1.6&range_c) = 1.6;
% %                 estEpsInf(estEpsInf>70&range_c) = 70; % for water
% %                 estCond(estCond> 1.59&range_c) = 1.59; % for water
%                  estEpsDelta(estEpsDelta~= 0&range_c) = 0;
                 
%for breast model
% values1=[4.68,23];
% values2=[3.12,33];
% values3=[0.0697,0.8];
% [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);       

values1=[1,100];
values2=[0,100];
values3=[0,100];
[estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3); 

    
        %============= plot all result in multiple figures===============
        if plot_flag
            Plot_All_Result;
        end
%         if iter_inner==15
%         figure()
%         imagesc(squeeze(imrotate(estEpsInf,90)));
%         caxis('auto')
%         colorbar;
%         axis image;
%         figure()
%         imagesc(squeeze(imrotate(estEpsDelta,90)));
%         caxis('auto')
%         colorbar;
%         axis image;
%         figure()
%         imagesc(squeeze(imrotate(estCond,90)));
%         caxis('auto')
%         colorbar;
%         axis image;
%         end
        %============= save the constructed parameters=================
        Result.L_curve_x(iter_total,1) = norm(contrast);
        Result.L_curve_Ax_b(iter_total,1) = norm(Mat_A * contrast - Data_ant);
        %   objective=rho; % this line is for CGLS
        pre_objective = log10(max(objective));
        Result.objective_all{iter_total} = objective;
        %     tau_all(iter_total)=tau_max;
        pre_contrast = contrast;
        iter_total = iter_total+1;
        
         Result.Epsinf{iter_inner,1} = estEpsInf;
    Result.EpsDelta{iter_inner,1} = estEpsDelta;
    Result.Cond{iter_inner,1} = estCond;
    end
    Result.hopping_result{ind_freq,1} = estEpsInf;
    Result.hopping_result{ind_freq,2} = estEpsDelta;
    Result.hopping_result{ind_freq,3} = estCond;
    
   
    
end % DBIM iterations
delete(FDTD_Prep_path);
end % end function