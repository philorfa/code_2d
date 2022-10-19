% LOAD VARIABLES FROM MIAN PROGRAM
% input: delT,pml_coarse,dimX,dimY,bbox_interior_mask,mu0£¬eps0£¬
% delX£¬db_immersion, antLocations_coarse_2D,
function  FDTD_Prep_path = FDTD_Preprocess(filename,delT, pml_coarse, dimX, dimY, bbox_interior_mask,...
    mu0, eps0, delX, EpsS, antLocations_coarse_2D)
% nmax =timeSteps;  
dt=delT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PML thickness in each direction
pml_size=pml_coarse;    % pml_coarse is loaded in sim_params_2D.mat
Imax = dimX+2*pml_size; 
Jmax = dimY+2*pml_size;
XBB=(1+pml_size):(Imax-pml_size);
YBB=(1+pml_size):(Jmax-pml_size);
% transfer "bbox_interior_mask" in main field to an new index for the extended field (include PML)
[I,J]=ind2sub([dimX dimY],bbox_interior_mask);
bbox_interior_mask_extend=sub2ind([Imax Jmax],I+pml_size,J+pml_size);
% transfer the location of the antenna to an extended field with PML
isource=antLocations_coarse_2D+pml_size; %the position of antennas
isource_ind=sub2ind([Imax, Jmax],isource(:,1),isource(:,2));
%======================set PML parameters===========================
eta = sqrt(mu0/(eps0*EpsS));
poly_m = 3;  % polynomial order for pml grading              
alpha_m=1.0;
m=poly_m;
ma=alpha_m;
sig_opt=0.8*(poly_m+1.0)/(eta*delX); 
sig_x_max = sig_opt;
sig_y_max =sig_x_max;
alpha_x_max = 0.03; %0.24;  
alpha_y_max = alpha_x_max; 
kappa_x_max = 1.0; %2.0;
kappa_y_max = kappa_x_max;
%===============================================
%================FDTD parameters for E and H fields====================
% E-Field
% x-dir
be_x_1=zeros(pml_size,1);
ce_x_1=zeros(pml_size,1);
alphae_x_PML_1=zeros(pml_size,1);
sige_x_PML_1=zeros(pml_size,1);
kappae_x_PML_1=zeros(pml_size,1);
be_x_2=zeros(pml_size,1);
ce_x_2=zeros(pml_size,1);
alphae_x_PML_2=zeros(pml_size,1);
sige_x_PML_2=zeros(pml_size,1);
kappae_x_PML_2=zeros(pml_size,1);
% y-dir
be_y_1=zeros(pml_size,1);
ce_y_1=zeros(pml_size,1);
alphae_y_PML_1=zeros(pml_size,1);
sige_y_PML_1=zeros(pml_size,1);
kappae_y_PML_1=zeros(pml_size,1);
be_y_2=zeros(pml_size,1);
ce_y_2=zeros(pml_size,1);
alphae_y_PML_2=zeros(pml_size,1);
sige_y_PML_2=zeros(pml_size,1);
kappae_y_PML_2=zeros(pml_size,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H-Field
% x-dir
bh_x_1=zeros(pml_size-1,1);
ch_x_1=zeros(pml_size-1,1);
alphah_x_PML_1=zeros(pml_size-1,1);
sigh_x_PML_1=zeros(pml_size-1,1);
kappah_x_PML_1=zeros(pml_size-1,1);
bh_x_2=zeros(pml_size-1,1);
ch_x_2=zeros(pml_size-1,1);
alphah_x_PML_2=zeros(pml_size-1,1);
sigh_x_PML_2=zeros(pml_size-1,1);
kappah_x_PML_2=zeros(pml_size-1,1);
% y-dir
bh_y_1=zeros(pml_size-1,1);
ch_y_1=zeros(pml_size-1,1);
alphah_y_PML_1=zeros(pml_size-1,1);
sigh_y_PML_1=zeros(pml_size-1,1);
kappah_y_PML_1=zeros(pml_size-1,1);
bh_y_2=zeros(pml_size-1,1);
ch_y_2=zeros(pml_size-1,1);
alphah_y_PML_2=zeros(pml_size-1,1);
sigh_y_PML_2=zeros(pml_size-1,1);
kappah_y_PML_2=zeros(pml_size-1,1);

% Denominators for the update equations 
den_ex=zeros(Imax,Jmax);
den_hx=zeros(Imax-1,Jmax-1);
den_ey=zeros(Imax,Jmax);
den_hy=zeros(Imax-1,Jmax-1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  SET CPML PARAMETERS IN EACH DIRECTION
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  x-dir
   for i = 1:pml_size
      sige_x_PML_1(i) = sig_x_max * ( (pml_size - i) / (pml_size - 1.0) )^m;
      alphae_x_PML_1(i) = alpha_x_max*((i-1.0)/(pml_size-1.0))^ma;
      kappae_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*((pml_size - i) / (pml_size - 1.0))^m;
      be_x_1(i) = exp(-(sige_x_PML_1(i) / kappae_x_PML_1(i) + alphae_x_PML_1(i))*dt/eps0);
      if ((sige_x_PML_1(i) == 0.0) && (alphae_x_PML_1(i) == 0.0) && (i == pml_size)) 
         ce_x_1(i) = 0.0;
      else
         ce_x_1(i) = sige_x_PML_1(i)*(be_x_1(i)-1.0)/(sige_x_PML_1(i)+kappae_x_PML_1(i)*alphae_x_PML_1(i))/ kappae_x_PML_1(i);
      end
   end
   for i = 1:pml_size-1
      sigh_x_PML_1(i) = sig_x_max * ( (pml_size - i - 0.5)/(pml_size-1.0))^m;
      alphah_x_PML_1(i) = alpha_x_max*((i-0.5)/(pml_size-1.0))^ma;
      kappah_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*((pml_size - i - 0.5) / (pml_size - 1.0))^m;
      bh_x_1(i) = exp(-(sigh_x_PML_1(i) / kappah_x_PML_1(i) + alphah_x_PML_1(i))*dt/eps0);
      ch_x_1(i) = sigh_x_PML_1(i)*(bh_x_1(i)-1.0)/(sigh_x_PML_1(i)+kappah_x_PML_1(i)*alphah_x_PML_1(i))/ kappah_x_PML_1(i);
   end
%-------------------
   for i = 1:pml_size
      sige_x_PML_2(i) = sig_x_max * ( (pml_size - i) / (pml_size - 1.0) )^m;
      alphae_x_PML_2(i) = alpha_x_max*((i-1.0)/(pml_size-1.0))^ma;
      kappae_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*((pml_size - i) / (pml_size - 1.0))^m;
      be_x_2(i) = exp(-(sige_x_PML_2(i) / kappae_x_PML_2(i) + alphae_x_PML_2(i))*dt/eps0);
      if ((sige_x_PML_2(i) == 0.0) && (alphae_x_PML_2(i) == 0.0) && (i == pml_size)) 
         ce_x_2(i) = 0.0;
      else
         ce_x_2(i) = sige_x_PML_2(i)*(be_x_2(i)-1.0)/(sige_x_PML_2(i)+kappae_x_PML_2(i)*alphae_x_PML_2(i))/ kappae_x_PML_2(i);
      end
   end
   for i = 1:pml_size-1
      sigh_x_PML_2(i) = sig_x_max * ( (pml_size - i - 0.5)/(pml_size-1.0))^m;
      alphah_x_PML_2(i) = alpha_x_max*((i-0.5)/(pml_size-1.0))^ma;
      kappah_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*((pml_size - i - 0.5) / (pml_size - 1.0))^m;
      bh_x_2(i) = exp(-(sigh_x_PML_2(i) / kappah_x_PML_2(i) + alphah_x_PML_2(i))*dt/eps0);
      ch_x_2(i) = sigh_x_PML_2(i)*(bh_x_2(i)-1.0)/(sigh_x_PML_2(i)+kappah_x_PML_2(i)*alphah_x_PML_2(i))/ kappah_x_PML_2(i);
   end
%  y-dir
   for j = 1:pml_size
      sige_y_PML_1(j) = sig_y_max * ( (pml_size - j ) / (pml_size - 1.0) )^m;
      alphae_y_PML_1(j) = alpha_y_max*((j-1)/(pml_size-1.0))^ma;
      kappae_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*((pml_size - j) / (pml_size - 1.0))^m;
      be_y_1(j) = exp(-(sige_y_PML_1(j) / kappae_y_PML_1(j) + alphae_y_PML_1(j))*dt/eps0);
      if ((sige_y_PML_1(j) == 0.0) && (alphae_y_PML_1(j) == 0.0) && (j == pml_size))
         ce_y_1(j) = 0.0;
      else
         ce_y_1(j) = sige_y_PML_1(j)*(be_y_1(j)-1.0)/(sige_y_PML_1(j)+kappae_y_PML_1(j)*alphae_y_PML_1(j))/ kappae_y_PML_1(j);
      end
   end
   for j = 1:pml_size-1
      sigh_y_PML_1(j) = sig_y_max * ( (pml_size - j - 0.5)/(pml_size-1.0))^m;
      alphah_y_PML_1(j) = alpha_y_max*((j-0.5)/(pml_size-1.0))^ma;
      kappah_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*((pml_size - j - 0.5) / (pml_size - 1.0))^m;
      bh_y_1(j) = exp(-(sigh_y_PML_1(j) / kappah_y_PML_1(j) + alphah_y_PML_1(j))*dt/eps0);
      ch_y_1(j) = sigh_y_PML_1(j)*(bh_y_1(j)-1.0)/(sigh_y_PML_1(j)+kappah_y_PML_1(j)*alphah_y_PML_1(j))/ kappah_y_PML_1(j);
   end
%-------------------
   for j = 1:pml_size
      sige_y_PML_2(j) = sig_y_max * ( (pml_size - j ) / (pml_size - 1.0) )^m;
      alphae_y_PML_2(j) = alpha_y_max*((j-1)/(pml_size-1.0))^ma;
      kappae_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*((pml_size - j) / (pml_size - 1.0))^m;
      be_y_2(j) = exp(-(sige_y_PML_2(j) / kappae_y_PML_2(j) + alphae_y_PML_2(j))*dt/eps0);
      if ((sige_y_PML_2(j) == 0.0) && (alphae_y_PML_2(j) == 0.0) && (j == pml_size))
         ce_y_2(j) = 0.0;
      else
         ce_y_2(j) = sige_y_PML_2(j)*(be_y_2(j)-1.0)/(sige_y_PML_2(j)+kappae_y_PML_2(j)*alphae_y_PML_2(j))/ kappae_y_PML_2(j);
      end
   end
   for j = 1:pml_size-1
      sigh_y_PML_2(j) = sig_y_max * ( (pml_size - j - 0.5)/(pml_size-1.0))^m;
      alphah_y_PML_2(j) = alpha_y_max*((j-0.5)/(pml_size-1.0))^ma;
      kappah_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*((pml_size - j - 0.5) / (pml_size - 1.0))^m;
      bh_y_2(j) = exp(-(sigh_y_PML_2(j) / kappah_y_PML_2(j) + alphah_y_PML_2(j))*dt/eps0);
      ch_y_2(j) = sigh_y_PML_2(j)*(bh_y_2(j)-1.0)/(sigh_y_PML_2(j)+kappah_y_PML_2(j)*alphah_y_PML_2(j))/ kappah_y_PML_2(j);
   end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  FILL IN DENOMINATORS FOR FIELD UPDATES
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ii = pml_size-1;
   for i = 1:Imax-1
      if (i <= pml_size-1)
         den_hx(i,:) = 1.0/(kappah_x_PML_1(i)*delX);
      elseif (i >= Imax+1-pml_size)
         den_hx(i,:) = 1.0/(kappah_x_PML_2(ii)*delX);
         ii = ii-1;
      else
         den_hx(i,:) = 1.0/delX;
      end
   end
   jj = pml_size-1;
   for j = 1:Jmax-1
      if (j <= pml_size-1) 
         den_hy(:,j) = 1.0/(kappah_y_PML_1(j)*delX);
      elseif (j >= Jmax+1-pml_size) 
         den_hy(:,j) = 1.0/(kappah_y_PML_2(jj)*delX);
         jj = jj-1;
      else
         den_hy(:,j) = 1.0/delX;
      end
   end
 %-------------------  
   ii = pml_size;
   for i = 1:Imax-1
      if (i <= pml_size)
         den_ex(i,:) = 1.0/(kappae_x_PML_1(i)*delX);
      elseif (i >= Imax+1-pml_size)
         den_ex(i,:) = 1.0/(kappae_x_PML_2(ii)*delX);
         ii = ii-1;
      else
         den_ex(i,:) = 1.0/delX;
      end
   end
   jj = pml_size;
   for j = 1:Jmax-1
      if (j <= pml_size)
         den_ey(:,j) = 1.0/(kappae_y_PML_1(j)*delX);
      elseif (j >= Jmax+1-pml_size)
         den_ey(:,j) = 1.0/(kappae_y_PML_2(jj)*delX);
         jj = jj-1;
      else
         den_ey(:,j) = 1.0/delX;
      end
   end
%%%======================temp value for computing====================
den_ex=den_ex(2:Imax-1,2:Jmax-1);
den_ey=den_ey(2:Imax-1,2:Jmax-1);
%---------- optimize the code on 12/11/2016--------
be_x_all=repmat(cat(1,be_x_1(2:end), zeros(Imax-2*pml_size,1),flip(be_x_2(2:end))) ,1,Jmax-2);
be_y_all=repmat(cat(2,be_y_1(2:end)',zeros(1,Jmax-2*pml_size),flip(be_y_2(2:end)')),Imax-2,1);
ce_x_all=repmat(cat(1,ce_x_1(2:end), zeros(Imax-2*pml_size,1),flip(ce_x_2(2:end))) ,1,Jmax-2);
ce_y_all=repmat(cat(2,ce_y_1(2:end)',zeros(1,Jmax-2*pml_size),flip(ce_y_2(2:end)')),Imax-2,1);

bh_x_all=repmat(cat(1,bh_x_1, zeros(Imax-1-2*(pml_size-1),1),flip(bh_x_2)),1,Jmax-1);
bh_y_all=repmat(cat(2,bh_y_1',zeros(1,Jmax-1-2*(pml_size-1)),flip(bh_y_2')),Imax-1,1);
ch_x_all=repmat(cat(1,ch_x_1, zeros(Imax-1-2*(pml_size-1),1),flip(ch_x_2)),1,Jmax-1);
ch_y_all=repmat(cat(2,ch_y_1',zeros(1,Jmax-1-2*(pml_size-1)),flip(ch_y_2')),Imax-1,1);
FDTD_Prep_path=['FDTD_Prep' filename '.mat'];
save(FDTD_Prep_path,'be_x_all','be_y_all','ce_x_all','ce_y_all','bh_x_all',...
    'bh_y_all','ch_x_all','ch_y_all','den_ex','den_ey','den_hy','den_hx','dt'...
    ,'Imax','Jmax','XBB','YBB','bbox_interior_mask_extend','isource_ind');
% [struc] = Combin_to_struc('FDTD_Prep.mat');
end