%% MRS signal Denosiing  using SCSA Method with splitting: 
% This function denoised the real part of the MRS signal using and optimal
% value of the h the ensures the preservation of the metabolite area,
% determined by <Peaks_area>, while remove as maximuim as possible the noise from the
% region determined by <Noise>

%% ######################  PARAMETERS ######################################
% Input 
% PRR    : The tolerated Peak Reduction ratio % used to check the
%            smoothness creteria ~5
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

function [yscsa, h_op, fs_op, Nh]=SCSA_Denoising_Total_Area( y, gm , fs)
% close all
N=max(size(y));

%% Sear optimal h for denoising
fprintf('\n-->Searching for the optimal h. Please wait...')

% get h bounds
 ymin=min(y);
 ymax=max(y);
 t0=0:1/fs:(N-1)/fs;
 Intg_y=trapz(sqrt(y-min(y)),t0);
 
h=abs(Intg_y);
fs0=floor(abs(Intg_y)/(10*ymax))+1;

cnt=0;

Results=[];
x_vec=[]; h_vec=[];fs_vec=[];
RMSE_vec=[];Cost_vec=[];
 for fs=[fs0:1:fs0*1000 ]%N/2] fs=[2:25 30:4:40]%N/2]
    cnt=cnt+1; x_vec(cnt)=fs;                    Results.x=x_vec;x_name='fs^*';


%  for h=linspace(0,Intg_y,N)
%     cnt=cnt+1;x_vec(cnt)=h;                    Results.x=x_vec;x_name='h^*';


    [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs,h,gm);
    
    h_vec=[h_vec h];           Results.h=h_vec;
    fs_vec=[fs_vec fs];        Results.fs=fs_vec;
    RMSE_vec=[RMSE_vec mse(y-yscsa)];        Results.RMSE=RMSE_vec;
    %% Compute cost function metric for denoising 
%     Peak_ratio=trapz(yscsa(Peaks_area))/trapz(y(Peaks_area));Cost_vec=[Cost_vec Peak_ratio];
    Peak_ratio=trapz(yscsa )/trapz(y );    Cost_vec=[Cost_vec Peak_ratio];
    Results.Totl_area=Cost_vec;
    
     
%% Get the optimal value 
% Totl_area
Idx1=find(Cost_vec==min(Cost_vec));
Smoot_op=Cost_vec(Idx1);
h_Totl_area_op=h_vec(Idx1);
fs_Totl_area_op=fs_vec(Idx1);

% RMSE
Idx2=find(RMSE_vec==min(RMSE_vec));
RMSE_op=RMSE_vec(Idx2);
h_RMSE_op=h_vec(Idx2);
fs_RMSE_op=fs_vec(Idx2);
    
    %  close all;
    figure(1);
      plot(y); hold on ;plot(yscsa);  hold off
      title(strcat('Setronstruction using : h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh)))

    figure(2);
      ax =plotyy(Results.x,Results.Totl_area ,Results.x,Results.RMSE);
      title(strcat('SCSA reconstruction using : h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh)...
                    ,'  , [Totl_area= ',num2str(Smoot_op),', h^*=', num2str(h_Totl_area_op), ', f_s^*=', num2str(fs_Totl_area_op),']',...
                    ', [RMSE=',num2str(RMSE_op),', h^*=', num2str(h_RMSE_op), ', f_s^*=', num2str(fs_RMSE_op), ']'))
        ylabel(ax(1), 'Totl_area');
        ylabel(ax(2), 'RMSE');
      pause(0.3)


    %% Stop if RMSE descreasing
    diff_RMSE=diff(RMSE_vec);
    dot2_RMSE=diff(diff_RMSE);
       
    Flag_RMSE=min(dot2_RMSE);
    if Flag_RMSE<0
        fs_op=fs_RMSE_op;
        h_op=h_RMSE_op;
        break;
    end
    
    
 end
 
 fs_op=fs_Totl_area_op;
 h_op=h_Totl_area_op;
 
 %% Get the optimal procesed signal 
     [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs_op,h_op,gm);
     
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
 
