function [t,zt0,f_ppm, zf0]=generate_Lorntzian()
fs=1;
N = 1024;% number of points 
amp = [10]*[ 1 2 0.3];
tau = [7.6]*[ 1 1.5 1.7];
freq = [1800]*[  0.5 1 2];

% Generation of the spectrum signal
f_Hz = 1:fs/N:fs;
t=0:1/fs:(N-1)/fs;
lorentz=f_Hz*0;
l=1:N;

for k=1:3
    
    lorentz = lorentz + amp(k)*(tau(k))./(1-i*((2*pi/N).*(freq(k)-4*l)).*(tau(k))); 
    
end

zf0=lorentz;
%% Get the time domaine signal 
zt0=ifft(zf0);

[f_ppm,f_Hz]= ppmScaleFid4(zt0, length(zt0));

