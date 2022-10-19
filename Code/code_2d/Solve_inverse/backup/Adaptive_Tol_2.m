function [Tol]=Adaptive_Tol_2(pre_objective)
    if pre_objective>2
        Tol=10^(-1);
    else if pre_objective>1.5
             Tol=10^(-2);
        else 
            Tol=10^(-3);
        end
    end
    disp( [ 'Tol=' num2str(Tol) ]);
    disp('adaptive tol is modified');
end
