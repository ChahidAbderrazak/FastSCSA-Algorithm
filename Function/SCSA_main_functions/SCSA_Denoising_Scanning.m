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

function [yscsa, h_op, fs_op, Nh,Noise_area]=SCSA_Denoising_Scanning( y, gm , fs ,y0)
%% List of paramters
PRR=0.08;              %  The tolerated Peak Reduction ratio % used to check the smoothness creteria ~ 0.08 
STD_TH_Coef=0.75;       % Paramter that determines the noisy region. decrese it to get wider  region 



% close all
N=max(size(y));

%% Get the noise and metabolites locations
[Noise_area]=Get_Noisy_Area_STD(y,STD_TH_Coef);
[Peaks_area]=Get_Peak_Area_STD(y);
% figure; plot(y); hold on;vline(Noise_area,'r');hold on;vline(Peaks_area,'g');hold off;

%% Sear optimal h for denoising
fprintf('\n-->Searching for the optimal h. Please wait...')

% get h bounds
 ymin=min(y);
 ymax=max(y);
 t0=0:1/fs:(N-1)/fs;
 Intg_y=trapz(sqrt(y-min(y)),t0);
 
h=Intg_y;
cnt=0;
Noise_area=1:300;

Results=[];
x_vec=[]; h_vec=[];fs_vec=[];
PSNR_vec=[];MSE_vec=[];Cost_vec=[];
 for fs=[2:25 30:4:40]%N/2]
    cnt=cnt+1; x_vec(cnt)=fs;                    Results.x=x_vec;x_name='fs^*';


%  for h=linspace(0,Intg_y,N)
%     cnt=cnt+1;x_vec(cnt)=h;                    Results.x=x_vec;x_name='h^*';


    [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs,h,gm);
    
    h_vec=[h_vec h];           Results.h=h_vec;
    fs_vec=[fs_vec fs];        Results.fs=fs_vec;
    MSE_vec=[MSE_vec mse(y-yscsa)];        Results.MSE=MSE_vec;
    PSNR_vec=[ PSNR_vec psnr(y0,yscsa)];     Results.PSNR=PSNR_vec;
    %% Compute cost function metric for denoising 
%     Peak_ratio=trapz(yscsa(Peaks_area))/trapz(y(Peaks_area));Cost_vec=[Cost_vec Peak_ratio];
%     Peak_ratio=trapz(yscsa )/trapz(y );    Cost_vec=[Cost_vec Peak_ratio];
    Smoothness=std(  diff(yscsa(Noise_area)))/abs(mean(diff(yscsa))  );    Cost_vec=[Cost_vec Smoothness];

    Results.Smooth=Cost_vec;
    
%     Amp_error=norm(y,yscsa)
     
%% Get the optimal value 
% Smooth
Idx1=find(Cost_vec==min(Cost_vec));
Smoot_op=Cost_vec(Idx1);
h_smooth_op=h_vec(Idx1);
fs_smooth_op=fs_vec(Idx1);


% Smooth
Idx2=find(PSNR_vec==max(PSNR_vec));
PSNR_op=Cost_vec(Idx2);
h_psnr_op=h_vec(Idx2);
fs_psnr_op=fs_vec(Idx2);


    
    %  close all;
    figure(1);
      plot(y); hold on ;plot(yscsa);  hold off
      title(strcat('Setronstruction using : h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh)))

    figure(2);
      ax =plotyy(Results.x,Results.Smooth ,Results.x,Results.PSNR);
      title(strcat('SCSA reconstruction using : h=',num2str(h),'  \gamma=',num2str(gm),'  , f_s= ',num2str(fs),'  , N_h=',num2str(Nh)...
                    ,'  , [Smooth= ',num2str(Smoot_op),', h^*=', num2str(h_smooth_op), ', f_s^*=', num2str(fs_smooth_op),']',...
                    ', [PSNR=',num2str(PSNR_op),', h^*=', num2str(h_psnr_op), ', f_s^*=', num2str(fs_psnr_op), ']'))
        ylabel(ax(1), 'Smooth');
        ylabel(ax(2), 'PSNR');
      pause(0.3)


 end
 
 fs_op=fs_smooth_op;
 h_op=h_smooth_op;
 
 %% Get the optimal procesed signal 
     [yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs_op,h_op,gm);

%  
% %% Search for the optimal h 
% for k=1:ref_loop
%         
%     fprintf('\n Refine %d : [%.2f,%.2f] .',ref_loop-k+1,h_min,h_max );
%     h_list= linspace(h_min,h_max,floor(nb_loops/2^k));
%     [h_op,Min_Cost]=search_4_optimal_h_denoising(t, Peaks_area, Noise_area,PRR, y,fs,h_list,gm);
%     h_min=h_op*(1-eps);
%     h_max=h_op*(1+eps);
%     %% plot the results     
%     h_op_list(k)=h_op;
%     Cost_function(k)=Min_Cost;
%     
% end

% h_op=min(h_op_list(find(Cost_function==min(Cost_function))));

% yscsa = scsa(h);

% if Padding==1
% 
%     [h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_padding(y,fs,h_op,gm);
% else
% 
%     [h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h_op,gm);
% 
% end



fprintf('\n--> Signal denoising is completed h=%f, hfs=%f,!!\n\n',h_op,fs_op)
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
    


function [noise_areas,Peaks_area]=Get_Noisy_Area_STD(y0,STD_TH_Coef)
N=max(size(y0));yf=y0;%(1:floor(0.48*N));

W=floor(N/50);  %the sliding windows size
for k=1:N-W-1
   std_y(k)=std(y0(k:W+k)); 
end

mean_std_y=1.4*mean(std_y)*(0*std_y+1);
% close all;plot(y0,'LineWidth',1.2); hold on;plot(std_y,'LineWidth',1.2); hold on;plot(mean_std_y,'LineWidth',1.2); hold on; legend('yf','std_y','mean_std_y') 

noise_areas=find(mean(std_y)>STD_TH_Coef*std_y);

%% Get a small noise region
A_noise0=floor(N/10);
if max(size(noise_areas))>A_noise0
    noise_areas=noise_areas(1:A_noise0);
end
 peaks=diff(noise_areas);
 peaks=[peaks 100];
    
a=1; cnt=0;
for k=1:size(peaks,2)
     if peaks(k)~=1 
         
         if k-a>5
            cnt=cnt+1;
            N_area(cnt,:)=[a k k-a];
            a=k+2;
         else
             a=k+1;
         end
     end
end

Idx=find(N_area(:,end)==max(N_area(:,end)));
noise_areas=N_area(Idx,1):N_area(Idx,2);
Peaks_areas=1:N;
Peaks_area=setdiff(Peaks_areas,noise_areas);
% close all; figure; plot(yf); hold on;  vline(noise_areas)
d=1;





function [Peaks_area]=Get_Peak_Area_STD(y)


%%Find Raw Peaks
idx=find(y==max(y));

Peaks_area=idx-3:idx+3;

% close all; figure; plot(y(peak_area)); hold on;  vline(noise_areas)
% 
% [pks,hrsAtpks] = findpeaks();
% plot(hrsAtpks,pks,'k^','MarkerSize',10)
% %%Find Peak Diff, (corresponding) Hr Diff
% pkDiff = abs(diff(pks));
% pkThresh = 50;
% pksTrans = pks(pkDiff>pkThresh);
% hrsTrans = hrsAtpks(pkDiff>pkThresh);
% plot(hrsTrans,pksTrans,'g^','MarkerSize',10)
% hold off
% legend('Original','Downsampled','Raw Peaks','Processed Peaks')
% ylim([-100 350])
% % % Debug PeakDiff if necessary for thresholding
% % plot(hrsAtpks(2:end),pkDiff,'r*')

% close all; figure; plot(yf); hold on;  vline(noise_areas)
d=1;




