% param=strcat(' h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh))
function plot_SCSA_denoising(figr, y, yscsa,param )
    %  close all;
    figure(figr);
        plot(y); hold on ;plot(yscsa);  hold off
        title(strcat('SCSA reconstruction [',param,']'))
        xlabel( 'Time');
        ylabel('Amplitude');

