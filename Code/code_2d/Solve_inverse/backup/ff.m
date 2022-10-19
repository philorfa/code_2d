function [x,out] = ff(A,b,alpha,opts)

%% Parameters and defaults
if isfield(opts,'lip'),    lip = opts.lip;     else lip = alpha*normest(A*A',1e-2); end
if isfield(opts,'tol'),    tol = opts.tol;     else tol = 1e-4;   end
if isfield(opts,'maxit'),  maxit = opts.maxit; else maxit = 500;  end
if isfield(opts,'maxT'),   maxT = opts.maxT;   else maxT = 1e3;   end
if isfield(opts,'x_ref'),  x_ref = opts.x_ref; else x_ref = []; out.out.hist_err = [];  end

n = size(A,2);

y = zeros(n,1);  % variable y in Nesterov's method
z = zeros(n,1);  % variable x in Nesterov's method
res = A'*b; % residual (b - Ax)
norm_b = norm(b);

shrink = @(z) sign(z).*max(0,abs(z)-alpha);

% extrapolation auxilary parameter, initialized to 1
t_pre=1;
t=(1+sqrt(1+4*t_pre^2))/2;  
%% Main iterations
start_time = tic;

for k = 1:maxit
    % --- extrapolation parameter ---
    beta = (t_pre-1)/t;   % computes beta_k in P.80 (2.2.9) of Nesterov's 2004 textbook
    
    % --- y-update ---
    z_new = y + res/lip;            % step 1a in P.80 of Nesterov's 2004 textbook
    y = z_new + beta*(z_new - z);   % step 1b in P.80 of Nesterov's 2004 textbook, extrapolation
    z = z_new;
    
    % --- x-update ---
    x = shrink(y);
    res = A'*b - A'*A*x; % res will be used in next y-update
    
    % --- theta-update ---
    t_pre=t;
    t = (1+sqrt(1+4*t_pre^2))/2;      % computes alp Accelerated Linearized Bregman Method, J. Sci. Comput, 2012. DOI: 10.1007/s10915-012-9592-9ha_{k+1} in P.80 (2.2.9) of Nesterov's 2004 textbook

    % --- diagnostics, reporting, stopping checks ---
    % reporting
    out.hist_obj(k) = norm(y) - norm(x); % dual objective
    out.hist_res(k) = norm(res); % residual size |b - Ax|_2
    if ~isempty(x_ref); out.hist_err(k) = norm(x - x_ref); end
    
    % stopping
    if toc(start_time) > maxT; break; end; % running time check
    if k > 1 && norm(res) < tol*norm_b; break; end % relative primal residual check
end

out.iter = k;

%% <http://www.caam.rice.edu/~optimization/linearized_bregman/accel/lbreg_accelerated.m Download this m-file>