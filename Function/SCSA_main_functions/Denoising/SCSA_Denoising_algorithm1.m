

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
% yf       : Noisy real art of the  complex MRS FID signal
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
function [y_denoised, h, Nh, split]=SCSA_Denoising_algorithm1(PRR, f_ppm, y_noisy, gm, fs)
global Plot_fig PRR  Noise_area    
Max_padd=0;
% [peak_pos, split]=find_peak_location(y_noisy);
N=max(size(y_noisy));split=[1 N 1]; peak_pos=1;

peak_pos=split(:,1);
figure; plot(y_noisy); hold on;vline(peak_pos,'r');hold on;%vline(Noise_area,'k');hold off ;

S_Sp=size(split,1);   
% Max_padd=200;

if Plot_fig==1
    figure(654); plot(y_noisy); hold on;
end

cntL=1;lgnd{cntL}='Noisy MRS signal';cntL=cntL+1;

for k=1:S_Sp
    
    if k==1
            Pad0=0;
            Pad1=Max_padd;
          
    elseif k==S_Sp
            Pad0=Max_padd;
            Pad1=0;
            
    else
            Pad0=Max_padd;
            Pad1=Max_padd;
                  
    end
    
    split_k=(split(k,1)-Pad0):(split(k,2)+Pad1);
    split_0=split(k,1):split(k,2);

    Coff=split(k,3);
    
    y=Coff*y_noisy(split_k);f_ppm_k=f_ppm(split_k);

    
    %% SCSA denoising using  PRR based algorithm 
    base=min(y);y=y-base;
    width_peaks=1;  % for split algorithm. this variable is not used 
    
    
%     [Noise_area,peaks_area]=find_Noisy_Area_STD(y);
    
%     Noise_area=Noisy_region(k,1):Noisy_region(k,2);
    
%     Noise_area=Noise_area-split_k(1)+1;
    
    [yscsa,h(k),Nh(k)]=SCSA_Denoising_PRR(PRR,f_ppm_k , y, gm , fs, width_peaks);
    yscsa=Coff*(yscsa(1+Pad0:end-Pad1) +base);
%     figure; plot(y); hold on;plot(yscsa);hold on;vline(Noise_area,'k');hold off ;
    
    y_denoised(split_0)=yscsa;
    
    %% plot sub_signals
    new_split_resultls(split_0)=yscsa;new_split_resultls(1:split_k(1)-1)=NaN;
    Noise_area_split=Noise_area+split_k(1);
    
    if Plot_fig==1
        figure(654);
        plot(new_split_resultls,'LineWidth',2); hold on ;    
        lgnd{cntL}=strcat('Denoised Sub-Signal',num2str(k));cntL=cntL+1;

        vline(Noise_area_split,'k');hold on ;
        lgnd{cntL}=strcat('A_',num2str(k));cntL=cntL+1;
    end

    clearvars new_split_resultls

    
%     A2=Noise_area_split;
end
  
if Plot_fig==1
    figure(654);
        A=legend(lgnd);
        A.FontSize=14;
        vline(peak_pos,'r');hold off ;
end
lgnd{cntL}=strcat('Denoised Sub-Signal',num2str(k));cntL=cntL+1;



%     figure; plot(y_noisy); hold on; hold on ;plot(y_denoised); hold on ;   

d=1;



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
% yf       : Noisy real art of the  complex MRS FID signal
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

function [yscsa, h_op, Nh,Noise_area]=SCSA_Denoising_PRR(PRR, f_ppm, yf, gm , fs , width_peaks)
% close all
global shif Plot_fig set_bound Bound_h rho Peaks_area
N=max(size(yf));

%% Get the noise and metabolites locations
% [Noise,Peaks_area]=find_metabolite_Noisy_Area(yf, width_peaks);
[Noise_area,Peaks_area]=find_Noisy_Area_STD(yf);
% figure; plot(yf); hold on;vline(Noise_area,'r');hold off;

% Peaks_areas=1:N;Noise_area=width_peaks;
% Peaks_area=setdiff(Peaks_areas,Noise_area);
% 
% Noise(2)=Noise(1);
%% Generate signals
%% MRS signal 
y = real(yf) ;%/max(real(yf))*76;

y_positive=y-min(y);

if set_bound==0
    h_max=(max(y_positive) + mean(y_positive));
    if max(y_positive)>1
        h_min=0.1;
    else
        h_min=max(y_positive)/10;
    end
    
else
    h_max=Bound_h(2);
    h_min=Bound_h(1);

end


h_max=h_max/rho;
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
            

function [noise_area,peaks_Locations]=find_metabolite_Noisy_Area(y0, width_peaks)
y0=y0-min(y0);
%% Find tw peaks location 
N=max(size(y0));yf=y0;%(1:floor(0.48*N));
TH=median(yf)+mean(yf(1:3));
yf_peaks=0*yf+TH;
yf_peaks(find(yf>TH))=yf(find(yf>TH));
yf_smooth   = sgolayfilt(yf_peaks, 2, 7);yf_smooth=yf_smooth-min(yf_smooth);
Th=max(yf_smooth) /4;
y_th=0*yf+Th;

idx=find(yf_peaks>Th);yf_peaks=0*yf;yf_peaks(idx)=yf(idx);

% close all; figure; plot(yf); hold on; plot(yf_peaks,'LineWidth',1.2); hold on;plot(y_th,'LineWidth',1.2); hold on; legend('yf','yf_peaks','y_th') 


[pks,locs,w,p] = findpeaks(yf_peaks);

peaks_location=locs(find(yf_peaks(locs) > 0.7*mean(yf_peaks(locs))));
Peaks_areas=0;

Peak_TH=peaks_location(find( yf_smooth(peaks_location)>Th));

for i=1:max(size(Peak_TH))
   
    peaks_Locations(i,:)=[Peak_TH(i)-width_peaks,Peak_TH(i)+width_peaks];    
end

if find(peaks_Locations<=0)>0
    
    peaks_Locations=peaks_Locations(2:end,:);
end


for i=1:size(peaks_Locations,1)
   
    Peaks_areas=[Peaks_areas,peaks_Locations(i,1):peaks_Locations(i,2)];
    
end

Peaks_area=1:peaks_location(end);

Peaks_areas=unique(Peaks_areas(2:end));
noise_areas=setdiff(Peaks_area,Peaks_areas);
noise_area=[noise_areas(1),noise_areas(end)];

% close all; figure; plot(yf(noise_area)); hold on; plot(yf_peaks,'LineWidth',1.2); hold on;plot(y_th,'LineWidth',1.2); hold on; legend('yf','yf_peaks','y_th') 

d=1;

function [noise_areas,Peaks_area]=find_Noisy_Area_STD0(y0, W)
N=max(size(y0));yf=y0;%(1:floor(0.48*N));

for k=1:N-W
   std_y(k)=std(y0(k:W+k)); 
end

mean_std_y=1.4*mean(std_y)*(0*std_y+1);
% close all; figure; plot(yf(noise_area)); hold on; 
% close all;plot(y0,'LineWidth',1.2); hold on;plot(std_y,'LineWidth',1.2); hold on;plot(mean_std_y,'LineWidth',1.2); hold on; legend('yf','std_y','mean_std_y') 

noise_areas=find(mean(std_y)>1.4*std_y);
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

% figure;plot(yf_peaks);hold on;plot(yf_smooth);hold on;plot(diff(yf_smooth));hold on;plot(yf_smooth*0+Th);


