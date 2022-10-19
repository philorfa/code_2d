function [Hx, Hy, psi_Hx, psi_Hy]=FDTD_Core_H(Ez,Hx, Hy, psi_Hx, psi_Hy,den_hy,den_hx,...
    DA,DB,bh_x_all,ch_x_all,bh_y_all,ch_y_all,delX)

temp_Ezx=-1*diff(Ez(1:end-1,:),1,2);
temp_Ezy=diff(Ez(:,1:end-1),1,1);
%  previous improved code 1
% % temp_Ezx=Ez-circshift(Ez,[0 -1]);
% % temp_Ezy=circshift(Ez,[-1 0])-Ez;
% % temp_Ezx=temp_Ezx(t_iH,t_jH);
% % temp_Ezy=temp_Ezy(t_iH,t_jH);

%  previous improved code 2
% % t_i=1:Imax-1;
% % t_j=1:Jmax-1;
% % temp_Ezx=Ez(t_i,t_j) - Ez(t_i,t_j+1);
% % temp_Ezy=Ez(t_i+1,t_j) - Ez(t_i,t_j);

Hx = DA * Hx + DB * (temp_Ezx .* den_hy);
psi_Hx = bh_y_all .* psi_Hx + ch_y_all .* temp_Ezx / delX;
Hx=Hx + DB * psi_Hx;

Hy = DA * Hy + DB *( temp_Ezy .* den_hx);
psi_Hy = bh_x_all .* psi_Hy + ch_x_all .* temp_Ezy / delX;
Hy = Hy + DB * psi_Hy;

end