function res_diff_fine=plot_diff_fine(Estimage,test_name,mul,res_diff_fine,iter_total)
load Debye_fine_062204_90slice_modified.mat
temp_image=imresize(Estimage,mul,'box');
res_diff_fine(iter_total) = norm(temp_image-Debye_fine_test_model_EpsInf)./...
    norm(Debye_fine_test_model_EpsInf);
figure(97)
plot(res_diff_fine(1:iter_total).','bo:');
title([ test_name '---difference from the fine']);
hold on, grid on, box on
end