%  %% Generate academic signal   
% [y0 t]=generate_signal();

%% Generate Lorntzian signal   
 [t,zt0,freq,zf0 ]=generate_Lorntzian();
 y0=real(zf0);

% Add Gaussian Noise  
y=awgn(y0,noise_level,'measured');



figure; plot(y0); hold on;  plot(y); hold off