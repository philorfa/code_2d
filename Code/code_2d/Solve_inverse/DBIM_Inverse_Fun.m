function Result = DBIM_Inverse_Fun(InP,total_freqs,maxIter,DBIM_mode,init_index,plot_flag,material_flag)
%  This is a matlab version of code for 2-D inverse scattering imaging
global multiple_sigma
% Original parameters
load sim_params_2D
% transformed model and relevant parameters.
[antLocations_coarse_2D, resolution, mul, coarse_grid_x, coarse_grid_y,...
    Debye_coarse_2D_model, numAnts] = Transform_Resolution(InP);%%%%%transform resolution from 0.5mm to selected (usually 2mm)
%=======disp number of ants=================
warning(['The number of antenna is ' num2str(numAnts)]);
%=====================initialize parameters=============================
timeSteps = 6000;% default timesteps. will be optimized by the function of Create_GMS
obs_iter = 5;  % number of iterations used to plot inverse solution
numTR = numAnts*(numAnts-1)/2;%%%%TR pair
PairTR = Generate_Ant_index(numAnts);

% computational domain dimensions
dimX = coarse_grid_x;%from temp_coarse_2D.mat
dimY = coarse_grid_y;%from temp_coarse_2D.mat
delX = resolution;%from temp_coarse_2D.mat
delT = delX/(2*c);
%  tauP=1.5e-11;
%=====================Create Debye model=============================
[estEpsInf, estEpsDelta, estCond, tauP, EpsS, bbSize, range_c,bbox_interior_mask, constraints]...
    = Load_FDTD_material_2D(Debye_coarse_2D_model, InP.model_phantom, ...
    DBIM_mode,material_flag);
%===============pre-process the FDTD parameters (CPML, etc.)===================
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
    [source,timeSteps]=create_GMS(fw.fctr,delT,'wide',timeSteps, dimX, dimY, 0);
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
        time_dbim_current=tic;
        % ============== run function for FDTD simulations ==============
        disp(['iteration #' num2str(iter_inner) ' > forward solver ...(FDTD is running)' ])
        %--------------------project_FDTD_2D-----------------------
        time_fdtd_current=tic;%%%fdtd time
        
        [greensFunctions_mag, greensFunctions_pha, receivedFields_mag,...
            receivedFields_pha, ~] = project_FDTD_2D(FDTD_Prep_path, estEpsInf, estEpsDelta, estCond,...
            tauP, delT, eps0, mu0, delX, source, numAnts, numFreqs, omegas,timeSteps, EpsS,...
            bbSize, esource_mag, esource_pha);%%%%%Note: greensFunctions_mag(pha) is not the real Green's fucntion.
        
        time_fdtd(iter_inner)=toc(time_fdtd_current);%%%fdtd time
        
        
        %=================inverse solution setup ==================================
        disp(['iteration #' num2str(iter_inner) ' > inverse solver ...'])
        [ Mat_A, Data_ant] = Generate_Matrix( receivedFields_mag, ...
            receivedFields_pha, greensFunctions_mag, greensFunctions_pha, measMag, ...
            measPha, numFreqs, numTR, omegas, eps0, PairTR, tauP, bbSize);
        %==========save residual and relative error before reconstruction=========
        % norm of residual
        for i = 1:numTR
            % store meas and calc antenna observations for TR pair
            mMag(i) = measMag(PairTR(i,1),PairTR(i,2),1);
            mPha(i) = measPha(PairTR(i,1),PairTR(i,2),1);
        end
        Result.residual3(iter_total) = norm(Data_ant,2);%/norm(mMag.*exp(1i.*mPha),2);
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
        %for experimental head model
        % values1=[20,100];
        % values2=[0,16.86];
        % values3=[0.1470,0.8];
        % [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        
        %for numerical breast model
        % values1=[2,24];
        % values2=[2,33];
        % values3=[0.0697,0.8];
        % [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        
        %for linear model
        % values1=[1,2];
        % values2=[0,0];
        % values3=[0.0055,0.055];
        % [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        
        %for head experiment
        %  values1=[20,78];
        %  values2=[15,100];
        %  values3=[0.1470,1.59];
        %  [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        
        %%for head experiment
        %  values1=[30,78];
        %  values2=[20,100];
        %  values3=[0.1470,1.59];
        %  [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        values1=[1,78];
        values2=[0,100];
        values3=[0.01,1.59];
        

        [estEpsInf,estEpsDelta,estCond]=SetConstraints(estEpsInf,estEpsDelta,estCond,range_c,values1,values2,values3);
        
        %         %============= plot all result in multiple figures===============
        if plot_flag
            Plot_All_Result;
        end
        %============= save the constructed parameters=================
        Result.L_curve_x(iter_total,1) = norm(contrast);
        Result.L_curve_Ax_b(iter_total,1) = norm(Mat_A * contrast - Data_ant);
        pre_objective = log10(max(objective));
        Result.objective_all{iter_total} = objective;
        iter_total = iter_total+1;
        
        Result.Epsinf{iter_inner,1} = estEpsInf;
        Result.EpsDelta{iter_inner,1} = estEpsDelta;
        Result.Cond{iter_inner,1} = estCond;
        time_dbim(iter_inner)=toc(time_dbim_current);
    end
    Result.hopping_result{ind_freq,1} = estEpsInf;
    Result.hopping_result{ind_freq,2} = estEpsDelta;
    Result.hopping_result{ind_freq,3} = estCond;
    Result.time_fdtd=time_fdtd;
    Result.time_dbim=time_dbim;
end % DBIM iterations
delete(FDTD_Prep_path);
end % end function