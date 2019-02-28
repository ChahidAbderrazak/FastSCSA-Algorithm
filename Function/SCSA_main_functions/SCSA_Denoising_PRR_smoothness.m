

%%  signal Denoising  using SCSA Method with splitting: 
% This function denoised the real part of the MRS signal using and optimal
% value of the h the ensures the preservation of the metabolite area,
% determined by <Peaks_area>, while remove as maximuim as possible the noise from the
% region determined by <Noise>

%% ######################  PARAMETERS ######################################
% Input 
% PRR    : The tolerated Peak Reduction ratio % used to check the
%            smoothness creteria ~5
% f_ppm   : The MRS spectrum fequancies in f_ppm
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

function [yscsa, h_op, Nh,Noise_area]=SCSA_Denoising_PRR_smoothness( f_ppm, y, gm , fs)


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

%% get h bounds
ymin=min(y)
ymax=max(y)
t0=0:1/fs:(N-1)/fs;
y_positive=y-min(y);
Intg_y=trapz(sqrt(y_positive),t0)
 
h_max=Intg_y;
h_min= max(y_positive)/100;

fprintf('\n-->Searching for the optimal h. Please wait...')

% you can choose either a vecotr for h or one value 
% h_list=h_min : .1 :h_max ;  % search in this interval
% h_list=10:0.1:17;
eps=0.2;
% shif=0;close all;


%% Setup the  Wait bar
global cnt_wait Tot_iter   wait_bar ref_loop nb_loops 
Tot_iter=0;
% nb_loops=  90;% 10%                % Number of loops 
% ref_loop =3;                       % number of refinement loop to find the optimal h for denoising
for k=1:ref_loop
    Tot_iter= Tot_iter+floor(nb_loops/2^k);
end

cnt_wait=0;
wait_bar = waitbar(cnt_wait,'Searching for the optimal value of h .... ','Name','MRS signal denoising using SCSA');


%% Search for the optimal h 
for k=1:ref_loop
        
    fprintf('\n Refine %d : [%.2f,%.2f] .',ref_loop-k+1,h_min,h_max );
    h_list= linspace(h_min,h_max,floor(nb_loops/2^k));
    [h_op,Min_Cost]=search_4_optimal_h_denoising(f_ppm, Peaks_area, Noise_area,PRR, y,fs,h_list,gm);
    h_min=h_op*(1-eps);
    h_max=h_op*(1+eps);
    %% plot the results     
    h_op_list(k)=h_op;
    Cost_function(k)=Min_Cost;
    
end

h_op=min(h_op_list(find(Cost_function==min(Cost_function))));

global    Noise_area  Padding
        % yscsa = scsa(h);
        
        if Padding==1

            [h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_padding(y,fs,h_op,gm);
        else

            [h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h_op,gm);
        
        end



fprintf('\n--> MRS denoising is completed h=%f!!\n\n',h_op)
close(wait_bar)


 %% Get the optimal     
% figure(2)
% plot(Cost_function); hold on
% plot(Cost_func_SNR); hold on
% plot(Cost_func_MSE); hold on
% 
% xlabel('Iterations')
% ylabel('Cost function ')
% title([ 'The Iterative MRS signal denoising']);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])


%% Plot the denoising results 
if Plot_fig==1
      figure;
    if shif==0
        plot(f_ppm,y,'b','LineWidth',2.5);hold on
    end
    
    plot(f_ppm, yscsa+5+shif ,'r','LineWidth',2);hold on
    plot(f_ppm, y-yscsa-5,'g','LineWidth',1.5); hold off
    
    legend({'Noisy input spectrum ', 'Denoised Spectrum ', 'Residual'},'Location','northwest');
    xlabel('ppm')
    ylabel('Intensity')
    set(gca,'YTickLabel',[])
    xlim([0 5])
    title([ ' SCSA MRS denoising with : h = ' num2str(h_op) , '   Nh = ' num2str(Nh) , ',   \beta = ' num2str(PRR) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
    set(gcf,'color','w') 
    set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
    text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
    box 
end
    


function  [h_op,Min_Cost]=search_4_optimal_h_denoising(f_ppm, Peaks_area, Noise_area,PRR, y,fs,h_list,gm)
global Tot_iter  wait_bar   cnt_wait  Plot_fig
h_list=h_list'; 
% start the loop for several values of h 
 for i=1:length(h_list)
 
    h = h_list(i);%0.045;%15.2;

    [Cost_SNR,Cost_peak,h, yscsa,Nh]=Cost_function_4_denoising(Peaks_area,Noise_area, y,fs,h,gm);
    Cost_func_SNR(i)=Cost_SNR;
    Cost_func_MSE(i)=Cost_peak;
    
    peaks=diff(Peaks_area);
    
    
    %% stop loop if you start going down
    idx00=find(y==max(y));
    idxxx=idx00-3:idx00+3;
    M(i)=norm(y(idxxx)-yscsa(idxxx));
    h_k(i)=h;
    
%     a=1; cnt=0;
%     PRR_peaks=0;
%     for k=1:size(peaks,2)
%      if    peaks(k)~=1 
%          
%          if k-a>5
%             M_area=a:k;
%             PRR_peaks=PRR_peaks+  abs(  max( abs(yscsa(M_area)))- max( abs(y(M_area)) ) );% / abs(max(abs(y(M_area))))
%             a=k+2;cnt=cnt+1;
%             
%          else
%              a=k+1;
%          end
%      end
%     end
%     
%     PRR_k(i,:)=PRR_peaks

    M_area=Peaks_area;
%     PRR_k(i,:)=( abs( max(yscsa(M_area))-max(y(M_area)) ) + abs( min(yscsa(M_area))-min(y(M_area)) ))/ abs(max(y(M_area)))
    PRR_k(i,:)=abs( max(yscsa(M_area))-max(y(M_area))  )/ abs(max(y(M_area)));
    
    if PRR_k(i)<PRR
    Cost_function(i)=Cost_SNR;%+Cost_peak;
    
    else
        
        Cost_function(i)=inf;
        
    end

    fprintf('.')
    cnt_wait=cnt_wait+1;
    % Update waitbar and message
    waitbar(cnt_wait/Tot_iter,wait_bar);
    
    if Plot_fig==1
    
        figure(154);
            plot(f_ppm,y,'LineWidth',2.5);hold on
            plot(f_ppm, yscsa ,'LineWidth',1.2);hold off

            % hold on
            % plot(f, y-yscsa-5,'g','LineWidth',1)
            legend({'Noisy input spectrum ', 'Denoised Spectrum '},'Location','northwest');
            xlabel('f_ppm')
            ylabel('Intensity')
            set(gca,'YTickLabel',[])
        %     xlim([0 5])
            title([ ' SCSA MRS densoing with : h = ' num2str(h) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
            set(gcf,'color','w') 
            set(gca,'Xdir','reverse');
            % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
            box 
            set(gcf, 'Position', get(0, 'Screensize'));            
            pause(0.3)
    end
    
 end
 
% Min_cost0=find(Cost_function==min(Cost_function));
% h_op=min(h_list(Min_cost0));



Cost_function=Cost_function';
Min_peak=find(M==min(M));

Min_Cost=min(Cost_function(Min_peak:end));
Min_cost_idx=find(Cost_function==Min_Cost);



h_op=min(h_list(Min_cost_idx));


%% plot the cost function 

% close all; figure;
%     plotyy(h_k, M, h_k, Cost_function);hold off
%         
%  d=1;

function [Cost_SNR,Cost_peak,h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]=Cost_function_4_denoising(Peaks_area,Noise_area, y,fs,h,gm)
global    Padding
        % yscsa = scsa(h);
        
        if Padding==1
            [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_padding(y,fs,h,gm);
        else
            
            [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);

        end


% coco1(i) = norm(y(a1:b1) - yscsa(a1:b1));
% coco2(i) =  norm(y(a2:b2) - yscsa(a2:b2));
% coco3(i) =  norm(y(a3:b3) - yscsa(a3:b3));
% % coco(i) =norm(y - yscsa);

y_res=y-yscsa;
Cost_peak=0;
for k=1:size(Peaks_area,1) 
    Cost_peak=Cost_peak + norm(y_res(Peaks_area(k,1:2)));
end


% SNR_y=max(yscsa)/std(y(Noise(1):Noise(2)));   %second term retaed to the  SNR
% SNR_scsa=max(yscsa)/std(yscsa(Noise(1):Noise(2)));   %second term retaed to the  residual SNR
% Cost_SNR=(1/abs(SNR_scsa));%*abs(SNR_y);

Cost_SNR=std(diff(yscsa(Noise_area)))/abs(mean(diff(yscsa)));
% Cost_SNR=1/Cost_SNR;
% Cost_SNR =std(y_res(Noise_area))   ;

% figure;plot(y_peaks);hold on;plot(y_smooth);hold on;plot(diff(y_smooth));hold on;plot(y_smooth*0+Th);


