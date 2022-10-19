function [estEpsInf, estEpsDelta, estCond, tauP, EpsS, bbSize, range_c,...
    bbox_interior_mask, constraints] = Load_FDTD_material_2D ...
    (Debye_coarse_2D_model, model_phantom, logic_skin,mode,init_index)
%bbox_interior_mask is the reconstruction area.
load(['..\data\model' num2str(model_phantom) '\all_material.mat']);
if (strcmp(mode,'Forward_original') || strcmp(mode,'Forward_background'))
    logic_skin=0;
end
if logic_skin % known skin in inverse process
    range_c=(Debye_coarse_2D_model.model==0.5);
    bbox_interior_mask=find(range_c);
else
    range_c=(Debye_coarse_2D_model.model==0.5 | Debye_coarse_2D_model.model==1|Debye_coarse_2D_model.model==4);
    bbox_interior_mask=find(range_c);
end
bbSize=length(bbox_interior_mask);


[dimX,dimY] = size(Debye_coarse_2D_model.model);
switch mode
    case 'Forward_original'
        estEpsInf   = Debye_coarse_2D_model.EpsInf;
        estEpsDelta = Debye_coarse_2D_model.DeltaEps;
        estCond     = Debye_coarse_2D_model.SigmaS;
        tauP        = mats.immer.fit.Tau;
        EpsS        = mats.immer.fit.EpsInf(1)+mats.immer.fit.DeltaEps(1);
        
    case {'Forward_background','Inverse'}
        % Note that the default of initial guess is the same as inner immersion.
        % Build 2D model based on the property of outer immersion
        estEpsInf   = mats.immer.fit.EpsInf*ones(dimX,dimY);
        estEpsDelta = mats.immer.fit.DeltaEps*ones(dimX,dimY);
        estCond     = mats.immer.fit.SigmaS*ones(dimX,dimY);
        tauP        = mats.immer.fit.Tau;
        EpsS        = mats.immer.fit.EpsInf+mats.immer.fit.DeltaEps;
        % insert the inner immersion

        %%%load materials
        [estEpsInf,estEpsDelta,estCond]=LoadMats(mats,Debye_coarse_2D_model.model,estEpsInf,estEpsDelta,estCond);
        
        %------------------------------------------------
        if strcmp(mode,'Inverse')
            estEpsInf(bbox_interior_mask)=init_guess.EpsInf(init_index);
            estEpsDelta(bbox_interior_mask)=init_guess.DeltaEps(init_index);
            estCond(bbox_interior_mask)=init_guess.SigmaS(init_index);
%             estEpsInf(bbox_interior_mask)=5.76; % match the original test setting
%             estEpsDelta(bbox_interior_mask)=5.51;% match the original test setting
%             estCond(bbox_interior_mask)=0.0802;% match the original test setting
        end
    otherwise
        error('The mode is not defined!');
end
% setup the property of tank on the model
if isfield(mats,'tank')
    estEpsInf   (Debye_coarse_2D_model.model==mats.tank.flag) = mats.tank.fit.EpsInf;
    estEpsDelta (Debye_coarse_2D_model.model==mats.tank.flag) = mats.tank.fit.DeltaEps;
    estCond     (Debye_coarse_2D_model.model==mats.tank.flag) = mats.tank.fit.SigmaS;
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