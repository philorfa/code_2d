function [E_fields_rec_dom_mag, E_fields_rec_dom_pha, receivedFields_mag,...
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
% eps_s = EpsS*ones(Imax, Jmax);
% eps_inf = EpsS*ones(Imax, Jmax);
% sigma_s = zeros(Imax, Jmax);
eps_s = zeros(Imax, Jmax);
eps_inf = zeros(Imax, Jmax);
sigma_s = zeros(Imax, Jmax);
% % % isource=antLocations_coarse_2D+pml_size; %the position of antennas
% % % isource_ind=sub2ind(size(sigma_s),isource(:,1),isource(:,2));
%------------------------------------------------------
eps_inf(XBB, YBB) = estEpsInf;
eps_s(XBB, YBB) = estEpsDelta + estEpsInf;
sigma_s(XBB, YBB) = estCond;
% background material
X_pml_min=min(XBB)-1;
X_pml_max=max(XBB)+1;
Y_pml_min=min(YBB)-1;
Y_pml_max=max(YBB)+1;

eps_inf(1:X_pml_min,1:Y_pml_min) = estEpsInf(1,1);
eps_inf(X_pml_max:end,Y_pml_max:end) = estEpsInf(end,end);
eps_inf(1:X_pml_min,Y_pml_max:end) = estEpsInf(1,end);
eps_inf(X_pml_max:end,1:Y_pml_min) = estEpsInf(end,1);
%%x-direction
eps_inf(XBB,1:Y_pml_min)=repmat( estEpsInf(:,1),1,Y_pml_min);
eps_inf(XBB,Y_pml_max:end)=repmat( estEpsInf(:,end),1,Jmax-Y_pml_max+1);
%%%y_direction
eps_inf(1:X_pml_min,YBB)=repmat( estEpsInf(1,:),X_pml_min,1);
eps_inf(X_pml_max:end,YBB)=repmat( estEpsInf(end,:),Imax-X_pml_max+1,1);

estEpsS=estEpsInf+estEpsDelta;
eps_s(1:X_pml_min,1:Y_pml_min) = estEpsS(1,1);
eps_s(X_pml_max:end,Y_pml_max:end) = estEpsS(end,end);
eps_s(1:X_pml_min,Y_pml_max:end) = estEpsS(1,end);
eps_s(X_pml_max:end,1:Y_pml_min) = estEpsS(end,1);
%%x-direction
eps_s(XBB,1:Y_pml_min)=repmat( estEpsS(:,1),1,Y_pml_min);
eps_s(XBB,Y_pml_max:end)=repmat( estEpsS(:,end),1,Jmax-Y_pml_max+1);
%%%y_direction
eps_s(1:X_pml_min,YBB)=repmat( estEpsS(1,:),X_pml_min,1);
eps_s(X_pml_max:end,YBB)=repmat( estEpsS(end,:),Imax-X_pml_max+1,1);


sigma_s(1:X_pml_min,1:Y_pml_min) = estCond(1,1);
sigma_s(X_pml_max:end,Y_pml_max:end) = estCond(end,end);
sigma_s(1:X_pml_min,Y_pml_max:end) = estCond(1,end);
sigma_s(X_pml_max:end,1:Y_pml_min) = estCond(end,1);
%%x-direction
sigma_s(XBB,1:Y_pml_min)=repmat( estCond(:,1),1,Y_pml_min);
sigma_s(XBB,Y_pml_max:end)=repmat( estCond(:,end),1,Jmax-Y_pml_max+1);
%%%y_direction
sigma_s(1:X_pml_min,YBB)=repmat( estCond(1,:),X_pml_min,1);
sigma_s(X_pml_max:end,YBB)=repmat( estCond(end,:),Imax-X_pml_max+1,1);

% eps_inf(eps_inf==0) = estEpsInf(1,1);
% eps_s(eps_s==0) = estEpsDelta(1,1) + estEpsInf(1,1);
% sigma_s(sigma_s==0) = estCond(1,1);
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
E_fields_rec_dom_mag = zeros(numAnts, bbSize, numFreqs);
E_fields_rec_dom_pha = zeros(numAnts, bbSize, numFreqs);
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
    temp_E_fields_imag = zeros(bbSize, numFreqs);
    temp_E_fields_real = zeros(bbSize, numFreqs);
    %====================BEGIN TIME STEP=======================
    for n = 1:timeSteps
        %=============update H====================
        [Hx, Hy, psi_Hx, psi_Hy] = FDTD_Core_H(Ez, Hx, Hy, psi_Hx, psi_Hy,...
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
        temp_E_fields_imag = temp_E_fields_imag + bboxObs_data * sin(omegas*n*dt);
        temp_E_fields_real = temp_E_fields_real + bboxObs_data * cos(omegas*n*dt);
%         if mod(n,100) == 0 
%       % ============== Plot current image ================ 
%         disp(n)
%         figure
%         temp = double(Ez(1:Imax,1:Jmax));
%         imagesc(squeeze(temp).');
%         cmax = max( max(max(temp)), -min(min(temp)) );
%         caxis([-cmax cmax])   
%         colorbar;
%         axis image;  
%         end
%       
    end %loop timesteps
    disp(['number of antenna is '  num2str(srcAnt)] )
    toc(scr);
    % convert observations to phasors and normalize by source
    temp_E_fields_imag = temp_E_fields_imag / (timeSteps/2);
    temp_E_fields_real = temp_E_fields_real / (timeSteps/2);
    % this is not the real green function, but is the Ez which will be used
    % to obtain green function later.
    E_fields_rec_dom_mag(srcAnt,:,:) = sqrt( temp_E_fields_imag.^2 + temp_E_fields_real.^2 ) ...
        ./ repmat(esource_mag,bbSize,1);
    E_fields_rec_dom_pha(srcAnt,:,:) = -atan2(temp_E_fields_imag,temp_E_fields_real) ...
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
