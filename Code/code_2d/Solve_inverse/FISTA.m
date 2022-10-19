function [y,objective,iternum] = FISTA(A,b,tau,varargin)
tolerance=1e-1;
maxIter=500;
minIter=1;
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch upper(varargin{i})
     case 'PSI'
       psi_function = varargin{i+1};
     case 'PHI'
       phi_function = varargin{i+1};   
     case 'STOPCRITERION'
       stopCriterion = varargin{i+1};
     case 'TOLERANCE'       
       tolerance = varargin{i+1};
     case 'MAXITER'
       maxIter = varargin{i+1};
     case 'MINITER'
       minIter = varargin{i+1};
     case 'INITIALIZATION'
       if prod(size(varargin{i+1})) > 1   % we have an initial x
	 init = 33333;    % some flag to be used below
	 x = varargin{i+1};
       else 
	 init = varargin{i+1};
       end
     case 'SPARSE'
       sparse = varargin{i+1};
     case 'TRUE_X'
       compute_mse = 1;
       true = varargin{i+1};
       if prod(double((size(true) == size(y))))
           plot_ISNR = 1;
       end
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end
  end
end
[~,n]=size(A);
x=zeros(n,1);
t=1;
L=1;
y=x;
AT= @(x) A'*x;
A = @(x) A*x;
phi_function = @(x) sum(abs(x(:)));
resid=A(y)-b;
pre_f=0.5*(resid'*resid) + tau*phi_function(y);
for iter=1:maxIter
    resid=A(y)-b;
    grad=AT(resid);
    stop_backtrack = 0 ;
    x_pre=x; 
%     k=0;
    while ~stop_backtrack
        x=soft(y-1/L*grad,tau/L);
        residx=b-A(x);
        temp1 = 0.5*(residx'*residx) ;
        temp2 = 0.5*(resid'*resid) + (x-y)'*grad + (L/2)*norm(x-y)^2 ;
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*2 ;
        end
%         k=k+1;
%         if k>1000
%             break;
%         end
    end
    f=temp1 + tau*phi_function(x);
%     if(norm(b-A*y)/norm(b)<tolerance)&&iter>5
%         condition=0;
%     end
    if(abs(f-pre_f)/pre_f<tolerance)&&(iter>minIter)
        break;
    end
    pre_f=f;
    objective(iter)=f;
    t_pre=t;
    t=(1+sqrt(1+4*t_pre^2))/2;
    y=x+(t_pre-1)/t*(x-x_pre);
end
iternum=iter
end

function y=soft(x,T)
y=sign(x).*max(0,abs(x)-T);
end