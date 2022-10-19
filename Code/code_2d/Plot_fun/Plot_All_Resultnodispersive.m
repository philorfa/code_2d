% Plot_All_Result
%%============norm of residual================================
figure(99)
plot(Result.residual3(1:iter_total).','ro:');
title([ InP.test_name '---residual for model' num2str(InP.model_phantom)]);
hold on, grid on, box on
refresh,pause(0.01)
%===============relative error========================
figure(98)
plot(Result.res_diff_1(1:iter_total).','bo:'); hold on
%plot(Result.res_diff_2(1:iter_total).','r*:'); hold on
plot(Result.res_diff_3(1:iter_total).','Gd:'); hold on
title([ InP.test_name '---relative error for model ' num2str(InP.model_phantom)]);
%legend('EpsInf','EpsDelta','Sigma_s');
legend('EpsInf','Sigma_s');
grid on, box on
refresh,pause(0.01)
%===============norm of difference from fine model==================
% %     res_diff_fine=plot_diff_fine(estEpsInf,test_name,mul,res_diff_fine,iter_total);
%==================== Plot current image ===========================
iter_plot=mod(iter_inner,obs_iter);
if iter_plot==0  %||iter_inner==1   estEpsInf
    H=ind_freq*10+ceil(iter_inner/20);
    sub_H=mod(iter_inner,20)/obs_iter;
    if sub_H==0
        sub_H=4;
    end
    Plot_DBIM(H,sub_H,estEpsInf,iter_inner,InP.test_name,fctrs)
end
%====================== add the original image===========================
if iter_inner==maxIter(ind_freq) && ind_freq==length(total_freqs) %
    if sub_H==4
        H=H+1;
        sub_H=1;
    else
        sub_H=sub_H+1;
    end
    Plot_DBIM(H,sub_H,Debye_coarse_2D_model.EpsInf,iter_inner,InP.test_name,fctrs)
end