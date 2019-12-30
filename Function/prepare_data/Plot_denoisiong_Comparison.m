

% param1= strcat('h^*=',num2str(h_op),' , f_s^*=',num2str(fs_op),', Nh= ',num2str(Nh_op) ,', PSNR=',num2str(PSNR_op))
function  Plot_denoisiong_Comparison(noise_level, f, yf0, yf,  yf_denoised1, param1, yf_denoised2, param2)

figure;
plot(f,yf,'k','LineWidth',2);hold on
plot(f,yf0,'g','LineWidth',1.7);hold on
plot(f,yf_denoised1,'r','LineWidth',1.4); hold on
plot(f,yf_denoised2,'b','LineWidth',1.3); hold on

legend({strcat('Noisy signal \sigma= ',num2str(noise_level)),strcat('Clean signal '),...
        strcat('Denoised signal', param1 ),strcat(' Denoised signal',param2 ) },...
        'Location','northeast');
% xlabel(' Hz')
% ylabel('Intensity')
% set(gca,'YTickLabel',[],'XTickLabel',[])
% xlim([f(1) f(end)])
title([ ' Pulse Shaped signal denoising using SCSA' ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
set(gcf,'color','w') 
% set(gca,'Xdir','reverse');
set(gca,'fontsize',16)

box 
