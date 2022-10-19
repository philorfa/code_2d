function [tau_TwIST]= Adaptive_tau(pre_contrast,percent)
        a=sort(abs(pre_contrast));
        b=round(length(pre_contrast)*percent);
        tau_TwIST=a(b);
        disp(['tau_TwIST=' num2str(tau_TwIST)])
end