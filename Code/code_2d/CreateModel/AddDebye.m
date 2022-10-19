function x=AddDebye(x,values)
    x.EpsInf=values(1);
    x.DeltaEps=values(2);
    x.SigmaS=values(3);
    if values(4)~=0
        x.Tau=values(4);
    end
end