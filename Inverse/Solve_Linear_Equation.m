%  =================== Run IMAT_CS for the inversion ====================
Anorm=norm(Mat_A,2);
Mat_A=Mat_A/Anorm;
Data_ant=Data_ant/Anorm;
%%max iteration
% iter=500;
% tol=1e-3;
%------------- test the adjustment factor
if numFreqs~=1
    Mat_A(:,2*end/3+1:end)=Mat_A(:,2*end/3+1:end)*multiple_sigma;
end
%-------------------------------------
switch InP.linear_method
    case 1
        %---------------------TwIST-----------------------
        tau_max=max(Mat_A'*Data_ant);
        if (ind_freq==1 && iter_inner==1)
            %tau_TwIST=0;
            Tol=10^(-1);
        else
            [Tol]=Adaptive_Tol_2(pre_objective);
        end
        %Tol=10^(-3);
        lam1=10^(-4);
        %Tol=tol;
        tau_TwIST=1e-2*tau_max;
        disp(['iteration #' num2str(iter_inner) ' >  TwIST...'])
        [contrast,~,objective,times,debias_start,mses,max_svd,iternum]=TwIST...
            (Data_ant,Mat_A,tau_TwIST,'TOLERANCEA',0,'lambda',lam1,'MAXITERA',100,'VERBOSE',0);
    case 2
        %--------------------------cgls algorithm-----------------------------
        disp(['iteration #' num2str(iter_inner) ' >  CGLS...'])
        [X,rho,eta] = cgls(Mat_A,Data_ant,25);
        contrast=X(:,25);
        objective=rho;
        iternum=25;
    case 3
        %--------------------------LSQR algorithm-----------------------------
        disp(['iteration #' num2str(iter_inner) ' >  LSQR...'])
        [contrast,FLAG,RELRES,iter] = lsqr(Mat_A,Data_ant,tol,iter);
        objective = RELRES;
    case 4  %%%%Direct solver
        disp(['iteration #' num2str(iter_inner) ' >  Tikhonov'])
        AA=Mat_A'*Mat_A;
        I=eye(size(AA));
        gamma=0.0001;
        b=Mat_A'*Data_ant;
        AAA=AA+gamma*I;
        contrast=AAA\b;
        objective = norm(b-AAA*contrast)/norm(b);
    case 5%%%%Nesterov
        disp(['iteration #' num2str(iter_inner) ' >  Nesterov'])
        contrast=nesterov(Mat_A,Data_ant,tol,iter);
        objective = norm(Data_ant-Mat_A*contrast)/norm(Data_ant);
    case 6 %%%FISTA
        tau_max=max(Mat_A'*Data_ant);
        %iter=50;
        iter=200;
        tol=1e-3;
        %disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
        [contrast,objective,iternum] = FISTA(Mat_A,Data_ant,1e-3*tau_max,'TOLERANCE',tol,'MAXITER',iter);
%     case 7
%         tau_max=max(Mat_A'*Data_ant);
%         iter=50;
%         disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
%         load('single_dbim.mat');
%         new_class=aaa;
%         new_class(new_class<1.5)=0.1;
%         new_class(new_class>=1.5)=1;
%         
%         for i=1:56
%             Mat_A(i,:)=Mat_A(i,:).*[new_class(bbox_interior_mask);new_class(bbox_interior_mask);new_class(bbox_interior_mask)]';
%         end
%         [contrast,objective,iternum] = FISTA(Mat_A,Data_ant,1e-2*tau_max,'TOLERANCE',0*tol,'MAXITER',iter);
%     case 8
%         tau_max=max(Mat_A'*Data_ant);
%         iter=50;
%         disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
%         load('single_both.mat');
%         new_class=bbb;
%         new_class(new_class==0)=0.2;
%         new_class(new_class==1)=0.5;
%         new_class(new_class==2)=1;
%         for i=1:56
%             Mat_A(i,:)=Mat_A(i,:).*[new_class(bbox_interior_mask);new_class(bbox_interior_mask);new_class(bbox_interior_mask)]';
%         end
%         [contrast,objective,iternum] = FISTA(Mat_A,Data_ant,1e-2*tau_max,'TOLERANCE',0*tol,'MAXITER',iter);
%     case 9
%         tau_max=max(Mat_A'*Data_ant);
%         iter=50;
%         disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
%         load('single_lsm.mat');
%         new_class=aaa;
%         new_class(aaa<1.5)=1e-1*tau_max;
%         new_class(aaa>=1.5)=1e-2*tau_max;
%         tau_fista=[new_class(bbox_interior_mask);new_class(bbox_interior_mask);new_class(bbox_interior_mask)]';
%         [contrast,objective,iternum] = FISTA1(Mat_A,Data_ant,tau_fista,'TOLERANCE',0*tol,'MAXITER',iter);
%     case 10
%         tau_max=max(Mat_A'*Data_ant);
%         iter=50;
%         disp(['iteration #' num2str(iter_inner) ' >  FISTA'])
%         load('single_both.mat');
%         new_class=bbb;
%         new_class(bbb==0)=1e-1*tau_max;
%         new_class(bbb==1)=1e-2*tau_max;
%         new_class(bbb==2)=1e-3*tau_max;
%         tau_fista=[new_class(bbox_interior_mask);new_class(bbox_interior_mask);new_class(bbox_interior_mask)]';
%         [contrast,objective,iternum] = FISTA1(Mat_A,Data_ant,tau_fista,'TOLERANCE',0*tol,'MAXITER',iter);
        
        
        
end
