%select the parameter of L1 regularization for TwIST
% This script includes three methods to choose an optimal parameter.
function L1_result=Optimization_L1(A,b,Tol,lam1,est_contrast,pre_solution,Original_Debye_model_EpsInf,visible)
%Original_Debye_model_EpsInf is set as Debye_coarse_test_model_EpsInf(bbox_interior_mask)
%set the range of value of regularization parameter
npoints = 200;  % Number of points on the L-curve for Tikh and dsvd.
smin_ratio = 16*eps;  % Smallest regularization parameter.
eta1 = zeros(npoints,1); rho1 = eta1; reg_param = eta1; RE1=eta1; Pt1=eta1;
reg_param2 = eta1; G=eta1;
%create regularization parameteer(method 1)
% % % [U,s,V]=csvd(A);
% % % reg_param(npoints) = max([s(end),s(1)*smin_ratio]);
% % % ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
% % % for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
%create regularization parameteer(method 2)
m=length(est_contrast);
sort_contrast=sort(abs(est_contrast));
for i=1:npoints
    point=ceil((i-1)/npoints*m);
    reg_param2(i)=sort_contrast(point+1);
end
reg_param2=flipud(reg_param2);
reg_param=reg_param2;
% First method(L1 curve)
% temp_ATA=trace(A'*A);
for i=1:npoints
    [contrast,x_debias,objective,times,debias_start,mses,max_svd]=TwIST...
          (b,A,reg_param(i),'TOLERANCEA',Tol,'lambda',lam1,'MAXITERA',5000);
     eta1(i)=norm(contrast,1);
     rho1(i)=norm(A*contrast-b);
     RE1(i)=norm((contrast(1:end/3)+pre_solution)-Original_Debye_model_EpsInf)./...
         norm(Original_Debye_model_EpsInf);
%      Pt1(i)=eta1(i)*rho1(i);
%      beta=(length(contrast)-temp_ATA)^2;
%      G(i)=rho1(i)^2/beta;
     solution{i}=contrast;
end
locate=find(RE1==min(RE1));
if visible
    figure(111); hold off;
    plot(eta1,rho1,'-*'); hold on;
    for i=20:20:200
        text(eta1(i),rho1(i),num2str(reg_param(i)));
    end
%     plot(eta1(locate),rho1(locate),'ro');
    figure(112); hold off;
    plot(reg_param,RE1,'-*'); hold on; xlabel('para'); ylabel('Relatvie error');
    figure(113); hold off;
    loglog(eta1,rho1,'-*'); hold on; xlabel('log(||x||_1)'); ylabel('log(||Ax-b||)');
    plot(eta1(locate),rho1(locate),'ro'); hold on;
    refresh,pause(0.1)
%     plot_lc(rho1,eta1,'-x',1,reg_param);
end
L1_result.para=reg_param;
L1_result.curve={eta1,rho1,RE1};
L1_result.solu=solution{locate};
L1_result.objective=0.5*(rho1(locate)).^2+reg_param(locate)*eta1(locate);
L1_result.RE_locate=locate;
end
