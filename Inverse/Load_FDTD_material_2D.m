function [estEpsInf, estEpsDelta, estCond, tauP, EpsS, bbSize, range_c,...
    bbox_interior_mask, constraints] = Load_FDTD_material_2D ...
    (Debye_coarse_2D_model, model_phantom, mode,material_flag)
%bbox_interior_mask is the reconstruction area.
load(['..\data\model' num2str(model_phantom) '\all_material.mat']);

%%%%%%reconstruction area range and size definition
step=size(material_flag);
range_c=(Debye_coarse_2D_model.model==material_flag(1));
for i=2:step(2)
    range_c=((Debye_coarse_2D_model.model==material_flag(i))|range_c);
end
bbox_interior_mask=find(range_c);
bbSize=length(bbox_interior_mask);
[dimX,dimY] = size(Debye_coarse_2D_model.model);
switch mode
    case 'Forward_original' %%%%%used to produce forward measurement data.
        estEpsInf   = Debye_coarse_2D_model.EpsInf;
        estEpsDelta = Debye_coarse_2D_model.DeltaEps;
        estCond     = Debye_coarse_2D_model.SigmaS;
        tauP        = mats.immer.fit.Tau;
        EpsS        = mats.immer.fit.EpsInf(1)+mats.immer.fit.DeltaEps(1);
        
    case {'Forward_background','Inverse'} %%%%used for reconstruction
        estEpsInf   = mats.immer.fit.EpsInf*ones(dimX,dimY);
        estEpsDelta = mats.immer.fit.DeltaEps*ones(dimX,dimY);
        estCond     = mats.immer.fit.SigmaS*ones(dimX,dimY);
        tauP        = mats.immer.fit.Tau;
        EpsS        = mats.immer.fit.EpsInf+mats.immer.fit.DeltaEps;
        %%%load materials
        [estEpsInf,estEpsDelta,estCond]=LoadMats(mats,Debye_coarse_2D_model.model,estEpsInf,estEpsDelta,estCond);
    otherwise
        error('The mode is not defined!');
end

disp(['The mode of loading material is ' mode]);

%set constraints there are some problems unsolved
if isfield(mats,'immer_inner')
    constraints.EpsInf = mats.immer_inner.fit.EpsInf;
    constraints.DeltaEps = mats.immer_inner.fit.DeltaEps;
    constraints.SigmaS = mats.immer_inner.fit.SigmaS;
else
    constraints.EpsInf = mats.immer.fit.EpsInf;
    constraints.DeltaEps = mats.immer.fit.DeltaEps;
    constraints.SigmaS = mats.immer.fit.SigmaS;
end
end