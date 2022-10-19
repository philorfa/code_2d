function [Ez_c,Jd,psi_Ezx,psi_Ezy] = FDTD_Core_Ez(Ez_c,Jd,psi_Ezx,psi_Ezy,CA,CB,Hy,Hx,den_ex,den_ey...
    ,Kd,Beta_d,dt,delX,be_x_all,ce_x_all,be_y_all,ce_y_all)
% t_i=2:Imax-1;
% t_j=2:Jmax-1;
% Ez_temp=Ez(t_i,t_j); % reduce the cost of indexing
Ez_former=Ez_c;
temp_Hy=diff(Hy(:,2:end),1,1);
temp_Hx=-1*diff(Hx(2:end,:),1,2);

%  previous improved code 1
% % temp_Hy=Hy-circshift(Hy,[1 0]);
% % temp_Hy=temp_Hy(t_iE,t_jE);

% % temp_Hx=circshift(Hx,[0 1])-Hx;
% % temp_Hx=temp_Hx(t_iE,t_jE);

% previous improved code 2
% % temp_Hy2=Hy(t_i,t_j) - Hy(t_i-1,t_j);
% % temp_Hx2=Hx(t_i,t_j-1) - Hx(t_i,t_j);
%=========================update main field of Ez==================================
Ez_c= CA.*Ez_c+CB.*(temp_Hy.*den_ex+temp_Hx.*den_ey - 0.5*(1+Kd)*Jd);
Jd=Kd*Jd+Beta_d.*(Ez_c-Ez_former)/dt;
%======================update pml of Ez in x direction====================
psi_Ezx = be_x_all.*psi_Ezx + ce_x_all.*temp_Hy/delX;
% Ez(t_i,t_j) = Ez(t_i,t_j) + CB(t_i,t_j).*psi_Ezx;
%======================update pml of Ez in y direction====================
psi_Ezy = be_y_all.*psi_Ezy + ce_y_all.*temp_Hx/delX;
% Ez(t_i,t_j) = Ez(t_i,t_j) + CB(t_i,t_j).*psi_Ezy;
Ez_c=Ez_c+CB.*(psi_Ezx+psi_Ezy);


%=========================update main field of Ez==================================
% % Ez_temp= CA.*Ez_temp+CB(t_i,t_j).*((Hy(t_i,t_j) - Hy(t_i-1,t_j)).*den_ex...
% %     +(Hx(t_i,t_j-1) - Hx(t_i,t_j)).*den_ey - 0.5*(1+Kd)*Jd);
% % Jd=Kd*Jd+Beta_d.*(Ez_temp-Ez_former)/dt;
% % % Ez(t_i,t_j)=Ez_temp;
% % %======================update pml of Ez in x direction====================
% % psi_Ezx = be_x_all.*psi_Ezx + ce_x_all.*(Hy(t_i,t_j) - Hy(t_i-1,t_j))/delX;
% % % Ez(t_i,t_j) = Ez(t_i,t_j) + CB(t_i,t_j).*psi_Ezx;
% % %======================update pml of Ez in y direction====================
% % psi_Ezy = be_y_all.*psi_Ezy + ce_y_all.*(Hx(t_i,t_j-1) - Hx(t_i,t_j))/delX;
% % % Ez(t_i,t_j) = Ez(t_i,t_j) + CB(t_i,t_j).*psi_Ezy;
% % 
% % Ez(t_i,t_j)=Ez_temp+CB(t_i,t_j).*(psi_Ezx+psi_Ezy);




end