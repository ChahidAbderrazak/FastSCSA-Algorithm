%% MRS signal Denosiing  using SCSA Method with splitting: 
% This function denoised the real part of the MRS signal using and optimal
% value of the h the ensures the preservation of the metabolite area,
% determined by <Peaks_area>, while remove as maximuim as possible the noise from the
% region determined by <Noise>

%% ######################  PARAMETERS ######################################
% Input 
% PRR    : The tolerated Peak Reduction ratio % used to check the
%            RMSEness creteria ~5
% t   : The MRS spectrum fequancies in t
% y       : Noisy real art of the  complex MRS FID signal
% Th_peaks_ratio  : The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
% width_peaks  : The with of the peak from its center(max values)
% gm, fs  : SCSA parameter: gm=0.5, fs=1

% Output
% yscsa: Densoied real art of the  complex MRS FID signal
% h_op : The  real art of the  complex MRS FID signal
% Nh   : Densoied real art of the  complex MRS FID signal

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: Oct,  2018
% King Abdullah University of Sciences and Technology (KAUST)

function [yscsa, h_op, fs_op, Nh]=SCSA_Denoising_PSNR( y, gm , fs ,y0)
% close all
N=max(size(y));

%% Sear optimal h for denoising
fprintf('\n-->Searching for the optimal h. Please wait...')

% get h bounds
 ymin=min(y);
 ymax=max(y);
 t0=0:1/fs:(N-1)/fs;
 Intg_y=abs(trapz(sqrt(y-min(y)),t0));
 
h=abs(Intg_y);

if h>100
    fs0=floor(h/ymax)+1;
else
    
    fs0=1;
    
end

Results=[];
x_vec=[]; h_vec=[];fs_vec=[];M_vec=[];
PSNR_vec=[];RMSE_vec=[];Cost_vec=[];MSE_vec=[];
cnt=0;

    
%%  for fs=[5:0.2:10]%N/2]
for fs=[fs0:4:fs0*1000 ]%N/2]
    cnt=cnt+1; x_vec(cnt)=fs;                    Results.x=x_vec;
    [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs,h,gm);M=(2*pi*h/fs);

%  for h=15:1:40
% %  for h=linspace(0,Intg_y,N)
%     cnt=cnt+1;x_vec(cnt)=h;                    Results.x=x_vec;
%     [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs,h,gm);M=(2*pi*h/fs);

%  for M=[0.01:0.05:2 ]%N/2]
%     cnt=cnt+1; x_vec(cnt)=M;                    Results.x=x_vec;
%     [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D_M(y,M,gm);

    M_vec=[M_vec M];           Results.M=M_vec;
    h_vec=[h_vec h];           Results.h=h_vec;
    fs_vec=[fs_vec fs];        Results.fs=fs_vec;
    RMSE_vec=[RMSE_vec mse(y-yscsa)];        Results.RMSE=RMSE_vec;
    PSNR_vec=[ PSNR_vec psnr(y0,yscsa)];     Results.PSNR=PSNR_vec;
    MSE_vec=[ MSE_vec mse(y0-yscsa)];     Results.MSE=MSE_vec;

%% Get the optimal value 
% RMSE
Idx1=find(RMSE_vec==min(RMSE_vec));
RMSE_op=RMSE_vec(Idx1);
h_RMSE_op=h_vec(Idx1);
fs_RMSE_op=fs_vec(Idx1);


% PSNR
Idx2=find(PSNR_vec==max(PSNR_vec));
PSNR_op=PSNR_vec(Idx2);
h_psnr_op=h_vec(Idx2);
fs_psnr_op=fs_vec(Idx2);


    %% Plot the denoised signal
    param=strcat(' h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh));
    figr=1;plot_SCSA_denoising(figr, y, yscsa,param )

    figure(20);
      ax =plotyy(Results.x,Results.RMSE ,Results.x,Results.PSNR);
      title(strcat('SCSA reconstruction using : h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh)...
                    ,'  , [RMSE= ',num2str(RMSE_op),', h^*=', num2str(h_RMSE_op), ', f_s^*=', num2str(fs_RMSE_op),']',...
                    ', [PSNR=',num2str(PSNR_op),', h^*=', num2str(h_psnr_op), ', f_s^*=', num2str(fs_psnr_op), ']'))
        ylabel(ax(1), '||y_{denoised} - y_{noisy}||');
        ylabel(ax(2), '||y_{denoised} - y_{clean}||');
        set(ax,'fontsize',16)

      pause(0.3)

    %% plot the first eigenfunction
    lgnd{cnt}=strcat('  f_s= ',num2str(fs));
    figure(21);%
    % plot(Df(500,:)); hold on
    subplot(411);plot(y0,'linewidth',2); 
    legend(lgnd)
    title('The input signal')
    
    subplot(412);plot(psinnor(:,1).^2,'linewidth',2); hold on
    legend(lgnd)
    title('\psi^2_1(f)')
    
    subplot(413);plot(psinnor(:,2).^2,'linewidth',2); hold on
    legend(lgnd)
    title('\psi^2_2(f)')

    subplot(414);plot(psinnor(:,3).^2,'linewidth',2); hold on
    legend(lgnd)
    title('\psi^2_3(f)')


      
      
      %% Stop if PSNR descreasing
      Flag_PSNR=min(diff(PSNR_vec));
      if Flag_PSNR<0
          
        fs_op=fs_psnr_op;
        h_op=h_psnr_op;

        break;
      end
      


 end


    
 %% Get the optimal procesed signal 
     [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs_op,h_op,gm);
     
%% Plot the denoised signal
    param=strcat(' h=',num2str(h_op),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs_op),'  , N_h=',num2str(Nh));
    figr=1;plot_SCSA_denoising(figr, y, yscsa,param )



fprintf('\n--> Signal denoising is completed h=%f, fs=%f,!!\n\n',h_op,fs_op)
% 
% 
% 
% %% Plot the denoising results 
% if Plot_fig==1
%       figure;
%     if shif==0
%         plot(t,y,'b','LineWidth',2.5);hold on
%     end
%     
%     plot(t, yscsa+5+shif ,'r','LineWidth',2);hold on
%     plot(t, y-yscsa-5,'g','LineWidth',1.5); hold off
%     
%     legend({'Noisy input spectrum ', 'Denoised Spectrum ', 'Residual'},'Location','northwest');
%     xlabel('ppm')
%     ylabel('Intensity')
%     set(gca,'YTickLabel',[])
%     xlim([0 5])
%     title([ ' SCSA MRS denoising with : h = ' num2str(h_op) , '   Nh = ' num2str(Nh) , ',   \beta = ' num2str(PRR) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
%     set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
%     set(gca,'fontsize',16)
%     text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
%     box 
% end
    
