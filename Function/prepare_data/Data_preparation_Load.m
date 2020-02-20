if sig_num==-3
    %% MRS signal 
    [f, y0_complex ]=generate_Lorntzian(noise_level); sig_name='Lorantzian_';
    yf0=real(y0_complex); N=max(size(yf0));

elseif sig_num==-2
    %% Load saved signal 
    load('./input/MRS_signal.mat') ;     yf0=y0;       sig_name='MRS_';

    
    
    elseif sig_num==-1
        %% Load saved signal 
        load('./input/pulse_signal.mat');    yf0=y0; f=t;   sig_name='Pulse_';  
            
else
    %% Academic signal
%      sig_num=6;   % 1)A1, 2) rect, 3)A1*t + A2*t, 4) A1*t^2 + A2*t^2, 5)A1*Exp(alpha*t)u + A2*Exp(beta*t)u, 
                 % 6)A1*sin(2*pi*alp.*t)+A2*sin(2*pi*beta.*t),  7) A1*sinc(alp.*t)+A2*sinc(beta.*t)
                 % 8) RC circuit decharge

    [yf0, f, name_sig]=generate_signal(sig_num);sig_name=strcat(name_sig,num2str(sig_num),'_');

                
                
end

%% Add Gaussian Noise 
% generate noise 
rng(1);N=max(size(yf0)); e= rand(size(yf0));e = e- mean(e);
yf = yf0 + noise_level*e; % noisy signal 


% %% plot
% close all;figure; plot(yf0); hold on;%  plot(yf); hold off

