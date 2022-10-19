function [source, Opt_timeSteps]=create_GMS(fc,delT,band_mode,timeSteps,dimX, dimY,test_mode)
% band=0.5e9;
% band=1.5e9;
disp(['Create a ' band_mode 'band source!'])
switch band_mode
    case 'wide'
        band=2e9;
    case 'narrow'
        band=1.0e9;
    otherwise
        error('source band is wrong');
end
fmin=fc-band;
fmax=fc+band;
omega=pi*(fmin+fmax);
Ts=2*2/(pi*(fmax-fmin));
% optimize Nc
Guass_limit=1e-5; %adjustment factor for finding an optimal delay length 
Nc=ceil(Ts/delT*sqrt(-1*log(Guass_limit)));

%Optimize the timesteps
Opt_timeSteps = ceil(2*Nc*1.5+max([dimX, dimY])*sqrt(2)*5);
if Opt_timeSteps>timeSteps
    Opt_timeSteps = timeSteps;
end
%Opt_timeSteps=1500;
disp(['The optimized timesteps is ' num2str(Opt_timeSteps)]);

N=1:Opt_timeSteps;
arg=(N-Nc)*delT/Ts;
source=exp(-1*arg.^2).*sin(omega*(N-1)*delT);
% source=sin(omega*N*delT);
% delay=5.7296e-10;
% tauS=1.9099e-10;
% source = sin(2*pi*fc.*((0:timeSteps-1)*delT-2*delay)) ...
%                      .* exp(-(((0:timeSteps-1)*delT-2*delay)./tauS).^2);
%plot(source)
if test_mode
    figure;
    plot(source);
    fft_image(source,1/delT);
end

end