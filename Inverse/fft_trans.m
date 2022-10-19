function [mag,pha]=fft_trans(source,omegas,delT,ref_mag,ref_pha)
[m,n]=size(source);
omegas=omegas(:)';

Re=source*cos((0:n-1)'*omegas*delT)/n*2;
Im=source*sin((0:n-1)'*omegas*delT)/n*2;
if ref_mag==0 & ref_pha==0
    mag=sqrt(Re.^2+Im.^2);
    pha=-atan2(Im,Re);
else
    temp_ref=repmat(ref_pha,m,1);
    mag=sqrt(Re.^2+Im.^2)*diag(1./ref_mag);
    pha=-atan2(Im,Re)-temp_ref;
end

end