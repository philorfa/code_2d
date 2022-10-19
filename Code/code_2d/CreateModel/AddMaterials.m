function mats=AddMaterials(mats,name,flag,values)
    fit=struct;
    fit=AddDebye(fit,values);
    f.fit=fit;
    f.flag=flag;
    mats.(name)=f;
end