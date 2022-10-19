%  =================== Run IMAT_CS for the inversion ====================
[m,n]=size(Mat_A);
B=zeros(m,n/3);
Mat_A=Mat_A(:,1:end/3);
Anorm=norm(Mat_A,2);
Mat_A=Mat_A/Anorm;
Data_ant=Data_ant/Anorm;
iter=5;%*mul;
%------------- test the adjustment factor
% if numFreqs~=1
%     Mat_A(:,2*end/3+1:end)=Mat_A(:,2*end/3+1:end)*multiple_sigma;
% end

%-------------------------------------
switch InP.linear_method
    case 1
        %---------------------TwIST-----------------------
        tau_max=max(Mat_A'*Data_ant);
        %tau_TwIST=0.1*max(Mat_A'*Data_ant);
        %[Tol]=Adaptive_Para(ind_freq);
        if (ind_freq==1 && iter_inner==1)
            %tau_TwIST=0;
            Tol=10^(-1);
        else
            [Tol]=Adaptive_Tol_2(pre_objective);
            %[tau_TwIST]=Adaptive_tau(pre_contrast,0.5);
        end   
        Tol=10^(-1);
        lam1=10^(-4);
        tau_TwIST=10^(-2)*tau_max;
        disp(['iteration #' num2str(iter_inner) ' >  TwIST...'])
        [contrast,~,objective,times,debias_start,mses,max_svd]=TwIST...
            (Data_ant,Mat_A,tau_TwIST,'TOLERANCEA',Tol,'lambda',lam1,'MAXITERA',50,'VERBOSE',0);
        if InP.Opt_L1
            L1_result(iter_total)=Optimization_L1(Mat_A,Data_ant,Tol,lam1,...
                contrast,estEpsInf(bbox_interior_mask),...
                Debye_coarse_2D_model.EpsInf(bbox_interior_mask),1);
            contrast=L1_result(iter_total).solu;
            objective=L1_result(iter_total).objective;
        end
    case 2
        %--------------------------cgls algorithm-----------------------------
        disp(['iteration #' num2str(iter_inner) ' >  CGLS...'])
        [X,rho,eta] = cgls(Mat_A,Data_ant,5);
        %[k_corner,info] = corner(rho,eta);
        %contrast=X(:,k_corner);
        %objective=k_corner;
         contrast=X(:,end);
         objective=rho;  
    case 3
        %--------------------------LSQR algorithm-----------------------------
        disp(['iteration #' num2str(iter_inner) ' >  LSQR...'])
        [contrast,FLAG,RELRES,iter] = lsqr(Mat_A,Data_ant,1e-1,50);
        objective = RELRES;
    case 4
        disp(['iteration #' num2str(iter_inner) ' >  Tikhonov'])
        AA=Mat_A'*Mat_A;  
        I=eye(size(AA));
        gamma=0.0001;
        b=Mat_A'*Data_ant;
        AAA=AA+gamma*I;
        contrast=AAA\b;
        objective = norm(b-AAA*contrast)/norm(b);
    case 5
        disp(['iteration #' num2str(iter_inner) ' >  Nesterov'])
        contrast=nesterov(Mat_A,Data_ant);
        objective = norm(Data_ant-Mat_A*contrast)/norm(Data_ant);
    case 6
        tau_max=max(Mat_A'*Data_ant);
        disp(['iteration #' num2str(iter_inner) ' >  Lbreg'])
        opts.maxit=iter;
        [contrast,out] = lbreg_accelerated(Mat_A,Data_ant,1e-2*tau_max,opts);
        objective = norm(Data_ant-Mat_A*contrast)/norm(Data_ant);
    case 7
        disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
         tau_max=max(Mat_A'*Data_ant);
        [contrast] = FISTA(Mat_A,Data_ant,1e-2*tau_max,1e-1,5);
        objective = norm(Data_ant-Mat_A*contrast)/norm(Data_ant);

end
Mat_A=[Mat_A,B,B];
contrast=[contrast;zeros(n/3,1);zeros(n/3,1)];