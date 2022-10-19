function [recondom,inf,delta,sigma]=ellipse(center,xr,yr,flag,recondom,inf,delta,sigma,setDebye,values)
[m,n]=size(inf);
for i=m/2-xr:m/2+xr
    for j=n/2-yr:n/2+yr
        if (i-center(1))^2/xr^2+(j-center(2))^2/yr^2<=1
            recondom(i,j)=flag;
            if setDebye==1
                inf(i,j)=values(1);
                delta(i,j)=values(2);
                sigma(i,j)=values(3);
            end
        end
    end
end

end