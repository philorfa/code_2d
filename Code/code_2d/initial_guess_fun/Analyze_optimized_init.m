function opt_init = Analyze_optimized_init(init_result,init_guess)
for i=1:length(init_result)
    temp_guess(i)=init_result(i).residual3(end); 
end
opt_init=find(temp_guess==min(temp_guess));
if length(opt_init)~=1
    warning('Existing more than one optimized initial guess!')
end
disp('the optimized initial guess is:');
disp(['EpsInf = ' num2str(init_guess.EpsInf(opt_init))]);
disp(['DeltaEps = ' num2str(init_guess.DeltaEps(opt_init))]);
disp(['EpsInf = ' num2str(init_guess.EpsInf(opt_init))]);
end