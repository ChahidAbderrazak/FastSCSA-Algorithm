function [frequency_ppm,ppmscale] = ppmScaleFid(u, nb_pts)
% Time vector (sampling)
dt= 8.3333e-04; % this is the 1200 Hz data from the 3T gent %8.66670e-004; sampling rate %.00024960; pour 9 tesla peut-etre 6.66670e-004; %
%8.3250
% 7.8290e-004 pour 3 Tesla
%dt = 0.0129;

%dt=1;

vt=[1:nb_pts].*dt;

% Creation of the axis in ppm
frequency=linspace(1/(2*dt),-1/(2*dt),nb_pts);   
% frequency=linspace(-1/(2*dt),1/(2*dt),nb_pts);
frequency_ppm=(frequency+600)/127.73;  %200.3; for 4.7 Tesla frequency 127.73 we added 600 which is equal to half of the BDW

% Matrice with the water resonance at 4.7 ppm instead of 0 ppm
u2=u'.*exp(-2i*pi*(4.7*127.73)*vt);   

%u2=u.*exp(-2i*pi*(-4.7*127.73).*vt);  conjugate data
%200.3 is the frequency of 4.7 Tesla
%127.7 is the frequency for 3 Tesla

% Fourier transform
%spectrum = real(fftshift(fft(u2(:,:))));
spectrum = abs((fft(conj(u2))));  %conj for wesi and csi data from Yao; it was real(fftshift(fft(conj(u2)))); 

frequency_ppm1 = circshift(frequency_ppm', -4);
spectrum1 = circshift(spectrum', 0);


% figure
% plot(frequency_ppm, spectrum);
% 
%  set(gca,'Xdir','reverse');
% % xlim([-5 +5])
%  grid on
ppmscale = spectrum;

