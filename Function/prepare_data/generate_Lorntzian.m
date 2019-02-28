function [f_ppm,zt0,zt,zf0,zf]=generate_Lorntzian(N)
amp = [10];
amp2 = [2.5]; 
tau = [7.6];
freq = [1800];
freq2 = [2500];


% the values ??of departure pp
%start_val = [1 16 1800 10 2500];
start_val = [20 10 1800 12 2500 ];

% Generation of the spectrum signal

N = 1024;% number of points 
T2 = 10;% time constant 
level = 5; % noise level 
aa = (1)^2;
f_Hz = 1:N;


for l=1: N 
    
    lorentz(l) = amp*(tau)/(1-i*((2*pi/N)*(freq-4*l))*(tau)) + amp2*(tau)/(1-i*((2*pi/N)*(freq2-4*l))*(tau)); 
    
end

zf0=lorentz;


% generate noise 
rng(1);
e= rand(1, 1024);
e = e- mean(e);

% generate the specturm 
zf = lorentz + level*e*(1+ i); % noisy signal 

zt=ifft(zf);
zt0=ifft(zf0);


[f_ppm,f_Hz]= ppmScaleFid4(zt, length(zt));

