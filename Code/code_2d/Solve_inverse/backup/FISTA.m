function [y] = FISTA(A,b,tau,tolerance,maxIter,varargin)
Atb=A'*b;
AtA=A'*A;
x=zeros(size(Atb));
x_pre=x;
iter=0;
t_pre=0;
t=1;
L=0.3*normest(A*A',1e-2);
y=x;
pre_f=0.5*norm(b-A*y)^2 + tau*(y'*y);
condition=1;
while iter<=maxIter &&condition
    iter=iter+1;
    y=x+(t_pre-1)/t*(x-x_pre);
    grad=AtA*y-Atb;
    stop_backtrack = 0 ;
    x_pre=x;
    f=0.5*norm(b-A*y)^2 + tau*(y'*y);
    while ~stop_backtrack
        x=soft(y-1/L*grad,tau/L);
        temp1 = 0.5*norm(b-A*x)^2 ;
        temp2 = 0.5*norm(b-A*y)^2 + (x-y)'*grad + (L/2)*norm(x-y)^2 ;
        if temp1 <= temp2
            stop_backtrack = 1 ;
        else
            L = L*2 ;
        end
    end
    if(norm(b-A*y)/norm(b)<tolerance)&&iter>5
        condition=0;
    end
    if(abs(f-pre_f)/pre_f<tolerance)&&iter>5
        condition=0;
    end
    pre_f=f;
    t_pre=t;
    t=(1+sqrt(1+4*t_pre^2))/2;
end


end

function y=soft(x,T)
y=sign(x).*max(0,abs(x)-T);
end