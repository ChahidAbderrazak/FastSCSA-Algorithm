% y0_complex: The signal without noise
% y_complex:The signal with noise

function [X,y0_complex,y_complex]=generate_Lorntzian(noise_level)
N = 1024;% number of points 
T2 = 10;% time constant 
amp = [10];
amp2 = [2.5];
tau = [7.6];
freq = [1800];
freq2 = [2500];

% the values ??of departure pp
%start_val = [1 16 1800 10 2500];
start_val = [20 10 1800 12 2500 ];

% Generation of the spectrum signal


%%%%%%%%%%%%
% noise_level = 10; % noise noise_level 
%%%%%%%%%%%

aa = (1)^2;
X = 1:N;


for l=1: N 
    
    lorentz(l) = amp*(tau)/(1-i*((2*pi/N)*(freq-4*l))*(tau)) + amp2*(tau)/(1-i*((2*pi/N)*(freq2-4*l))*(tau)); 
    
end

y0_complex=lorentz'; %% Th signal without noise 


%% generate noise 
rng(1);
e= rand(1, N);
e = e- mean(e);

%% generate the specturm 
y_complex = lorentz + noise_level*e+ noise_level*i*e; % noisy signal 
y_complex=y_complex';
