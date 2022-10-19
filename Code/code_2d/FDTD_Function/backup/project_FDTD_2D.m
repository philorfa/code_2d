function [greensFunctions_mag, greensFunctions_pha, receivedFields_mag,...
    receivedFields_pha, antObs_data_test] = project_FDTD_2D(FDTD_Prep_path, estEpsInf, ...
    estEpsDelta, estCond, tauP, dt, eps0, mu0, delX, source, numAnts, numFreqs, ...
    omegas, timeSteps, EpsS, bbSize, esource_mag, esource_pha)
%  This is a matlab version of code for 2-D Debye model FDTD with CPML
%  written by Zhenzhuang Miao (Supervisor: Dr. Panagiotis Kosmas)
%  04/09/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
running_time=tic;
load(FDTD_Prep_path);
%===============================================
eps_s = EpsS*ones(Imax, Jmax);
eps_inf = EpsS*ones(Imax, Jmax);
sigma_s = zeros(Imax, Jmax);
% % % isource=antLocations_coarse_2D+pml_size; %the position of antennas
% % % isource_ind=sub2ind(size(sigma_s),isource(:,1),isource(:,2));
%------------------------------------------------------
eps_inf(XBB, YBB) = estEpsInf;
eps_s(XBB, YBB) = estEpsDelta + estEpsInf;
sigma_s(XBB, YBB) = estCond;
% background material
eps_inf(eps_inf==0) = EpsS;
eps_s(eps_s==0) = EpsS;
sigma_s(sigma_s==0) = 0;
% add metal material for the point sources (important)
% sigma_s(sub2ind(size(eps_inf),isource(:,1),isource(:,2)))=3.27e7; %Aluminium
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  FILL IN UPDATING COEFFICIENTS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kd = (2 * tauP - dt) / (2 * tauP + dt);
Beta_d = (2 * eps0 * (eps_s - eps_inf) * dt) / (2 * tauP + dt);
DA = 1.0;
DB = (dt / mu0); % ????
CA = (2 * eps0 * eps_inf - sigma_s * dt + Beta_d) ./ ...
    (2 * eps0 * eps_inf + sigma_s * dt + Beta_d);
% % % compulsory the value at the antennas equal to 1, because that point
% % % source is not a Debye model.
CA = CA(2:Imax-1, 2:Jmax-1);
CB = 2*dt./(2*eps0.*eps_inf+sigma_s*dt+Beta_d);
CB = CB(2:Imax-1, 2:Jmax-1);
Beta_d = Beta_d(2:Imax-1, 2:Jmax-1); %temp value for computing

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%% PRE-ALLOCATE MEMORY (large chunks first) %%%%%
% for saving phasors inside FDTD observation bounding box
greensFunctions_mag = zeros(numAnts, bbSize, numFreqs);
greensFunctions_pha = zeros(numAnts, bbSize, numFreqs);
% for saving phasors at the antennas in FDTD
% % receivedFields_mag = zeros(numAnts,numAnts,numFreqs);
% % receivedFields_pha = zeros(numAnts,numAnts,numFreqs);
antObs_data_test = zeros(numAnts, numAnts, timeSteps);
% t_iH=1:Imax-1; %para for updating function H
% t_jH=1:Jmax-1; %para for updating function H
t_iE=2:Imax-1; %para for updating function E
t_jE=2:Jmax-1; %para for updating function E
%================= BEGIN TO SIMULATE ALL ANTENNAS====================
for srcAnt= 1:numAnts
    scr=tic;
    %=================initialize values for field======================
    %  BEGIN TIME STEP
    Ez = zeros(Imax, Jmax);
    Hx = zeros(Imax-1,Jmax-1);
    Hy = zeros(Imax-1,Jmax-1);
    Jd = zeros(Imax-2, Jmax-2);% t_i=2:Imax-1;    t_j=2:Jmax-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  CPML
    psi_Ezx = zeros(Imax-2, Jmax-2);
    psi_Ezy = zeros(Imax-2, Jmax-2);
    psi_Hx  = zeros(Imax-1, Jmax-1);
    psi_Hy  = zeros(Imax-1, Jmax-1);
    % temp storage for saving data inside observation bounding box
    tempGre_imag = zeros(bbSize, numFreqs);
    tempGre_real = zeros(bbSize, numFreqs);
    %====================BEGIN TIME STEP=======================
    for n = 1:timeSteps
        %=============update H====================
        [Hx, Hy, psi_Hx, psi_Hy] = FDTD_Core_H(Ez, Hx, Hy, psi_Hx, psi_Hy,...]
            den_hy, den_hx,DA,DB, bh_x_all, ch_x_all, bh_y_all, ch_y_all, delX);
        %=============UPDATE Hx====================
        % % %         [Hx,psi_Hx]=FDTD_Core_Hx(Ez,Hx,psi_Hx,Imax,Jmax,den_hy,DA,DB,bh_y_all,ch_y_all,delX);
        %=============UPDATE Hy====================
        % % %         [Hy,psi_Hy]=FDTD_Core_Hy(Ez,Hy,psi_Hy,Imax,Jmax,den_hx,DA,DB,bh_x_all,ch_x_all,delX);
        %=============UPDATE Ez====================
        [Ez(t_iE,t_jE),Jd,psi_Ezx,psi_Ezy] = FDTD_Core_Ez(Ez(t_iE, t_jE),Jd,...
            psi_Ezx, psi_Ezy, CA, CB, Hy, Hx, den_ex, den_ey, Kd, Beta_d,...
            dt, delX, be_x_all, ce_x_all, be_y_all, ce_y_all);
        % %         [Ez,Jd] = FDTD_Core_Ez_1(Ez,Jd,Imax,Jmax,CA,CB,Hy,den_ex,Hx,den_ey,Kd,Beta_d,dt);
        % %         [Ez,psi_Ezx,psi_Ezy]=FDTD_Core_Ez_2(Ez,psi_Ezx,psi_Ezy,Imax,Jmax,be_x_all,...
        % %     ce_x_all,Hy,Hx,delX,CB,be_y_all,ce_y_all);
        %====================source===========================
           %Ez(isource_ind(srcAnt)) = source(n);
              Ez(isource_ind(srcAnt)) = source(n) + Ez(isource_ind(srcAnt));
%         imagesc(Ez);
        % observe field data at every obstime_ds timesteps
        bboxObs_data=Ez( bbox_interior_mask_extend);
        antObs_data_test(srcAnt,:,n)=Ez(isource_ind);
        
        % Computing frequency-domain quantities by FFT
        tempGre_imag = tempGre_imag + bboxObs_data * sin(omegas*n*dt);
        tempGre_real = tempGre_real + bboxObs_data * cos(omegas*n*dt);
        %obstime=200;
        %if mod(n,obstime)==0
         %   disp(n)
          %  figure
           % temp = double(Ez(XBB,YBB));
            %Cmax = abs(max(max(max(max(temp))),min(min(min(temp)))));
            %color =[-Cmax Cmax];
            %imagesc(temp(:,:)');
            %title(Cmax);
            %caxis('auto')
            %caxis(color)
            %axis image;
        %end
    end %loop timesteps
    disp(['number of antenna is '  num2str(srcAnt)] )
    toc(scr);
    % convert observations to phasors and normalize by source
    tempGre_imag = tempGre_imag / (timeSteps/2);
    tempGre_real = tempGre_real / (timeSteps/2);
    % this is not the real green function, but is the Ez which will be used
    % to obtain green function later.
    greensFunctions_mag(srcAnt,:,:) = sqrt( tempGre_imag.^2 + tempGre_real.^2 ) ...
        ./ repmat(esource_mag,bbSize,1);
    greensFunctions_pha(srcAnt,:,:) = -atan2(tempGre_imag,tempGre_real) ...
        - repmat(esource_pha,bbSize,1);
end % loop over antenna simulations

[mag,pha] = fft_trans(reshape(antObs_data_test,numAnts^2,timeSteps),omegas,...
    dt,esource_mag,esource_pha);
receivedFields_mag = reshape(mag, numAnts, numAnts, numFreqs);
receivedFields_pha = reshape(pha, numAnts, numAnts, numFreqs);
disp('the whole running time of FDTD is ')
toc(running_time);
end
%***********************************************************************
%     END TIME-STEPPING LOOP
%**********************************************************************
