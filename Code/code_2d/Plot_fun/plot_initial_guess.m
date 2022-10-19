function plot_initial_guess(init_result,init_guess)
num_plot=length(init_result);

figure(100);
for i=1:num_plot
    plot(init_result(i).residual3);
    hold on
end
Add_Marker(100);
title('residual for initial guess');
leg=cat(2,repmat('EpsInf=',num_plot,1),num2str(init_guess.EpsInf(1:num_plot)'));
legend(leg);
grid on
xlabel('iteration');
ylabel('residual');


figure(101);
for i=1:num_plot
    plot(init_result(i).res_diff_1);
    hold on
end
Add_Marker(101);
title('EpsInf for initial guess');
leg=cat(2,repmat('EpsInf=',num_plot,1),num2str(init_guess.EpsInf(1:num_plot)'));
legend(leg);
grid on
xlabel('iteration');
ylabel('difference');

figure(102);
for i=1:num_plot
    plot(init_result(i).res_diff_2);
    hold on
end
Add_Marker(102);
title('EpsDelta for initial guess');
leg=cat(2,repmat('EpsDelta=',num_plot,1),num2str(init_guess.DeltaEps(1:num_plot)'));
legend(leg);
grid on
xlabel('iteration');
ylabel('difference');

figure(103);
for i=1:num_plot
    plot(init_result(i).res_diff_3);
    hold on
end
Add_Marker(103);
title('SigmaS for initial guess');
leg=cat(2,repmat('SigmaS=',num_plot,1),num2str(init_guess.SigmaS(1:num_plot)'));
legend(leg);
grid on
xlabel('iteration');
ylabel('difference');
end