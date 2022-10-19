function x=nesterov(A,b,tol,maxiter)
dim=size(A);
x=zeros(dim(2),1);
epsilon=0.01;
lamb_prev=0;
lamb_curr=1;
gamma=1;
y_prev=x;
L=normest(A*A',1e-2);
alpha=1/L;
gd=grad(A,b,x);
iter=0;
f_pre=norm(A*x-b,2);
while norm(gd)>=0.001*epsilon &&iter<=maxiter
    iter=iter+1;
    y_curr=x-alpha.*gd;
    x=(1+gamma)*y_curr-gamma*y_prev;
    y_prev=y_curr;
    lamb_temp=lamb_curr;
    lamb_curr=(1+sqrt(1+4*lamb_prev*lamb_prev))/2;
    lamb_prev=lamb_temp;
    gamma=(-1+lamb_curr)/lamb_prev;
    gd=grad(A,b,x);
    f=norm(A*x-b,2);
    if abs(f-f_pre)/f_pre<tol
        break;
    end
    f_pre=f;
end

end
function gd=grad(A,b,x)
gd=A'*(A*x-b);
end