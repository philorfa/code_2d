function [epsinf,epsdelta,sigma]=SetConstraints(epsinf,epsdelta,sigma,ind,values1,values2,values3)
epsinf(epsinf<values1(1)&ind)=values1(1);
epsinf(epsinf>values1(2)&ind)=values1(2);
epsdelta(epsdelta<values2(1)&ind)=values2(1);
epsdelta(epsdelta>values2(2)&ind)=values2(2);
sigma(sigma<values3(1)&ind)=values3(1);
sigma(sigma>values3(2)&ind)=values3(2);
end