function [Tol]=Adaptive_Para(ind_freq)
switch ind_freq
    case 1
        Tol=10^(-1);
    case 2
        Tol=10^(-1);
    case 3
        Tol=10^(-2);
    case 4
        Tol=10^(-2);
    case 5
        Tol=10^(-3);
    case 6
        Tol=10^(-3);
end
        

end