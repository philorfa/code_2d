function [estEpsInf,estEpsDelta,estCond]=LoadMats(mats,model,estEpsInf,estEpsDelta,estCond)
name=fieldnames(mats);
[m,~]=size(name);
for i=1:m
    disp(['add ' name{i} ' to the materials!'])
    estEpsInf(model==mats.(name{i}).flag) = mats.(name{i}).fit.EpsInf;
    estEpsDelta(model==mats.(name{i}).flag) = mats.(name{i}).fit.DeltaEps;
    estCond(model==mats.(name{i}).flag) = mats.(name{i}).fit.SigmaS;
end
end