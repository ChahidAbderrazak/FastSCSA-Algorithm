

%%
    %**********************************************************************
    %                           SCSA Densoing Function                            *
    %**********************************************************************

 
 %% Description
 % Script search for the h values that  densoie optimally the input noisy
 % signal y. the output signal is reconstructed using K-Nh eigenfunction
 % such that K-Nh= Nh0/ K. where Nh0 are the number of eigenfunctions needed
 % to reconstruct the noisy signal completelly.
 % Done: Feb,  2019

%% input parapeters
% K   :  The K-SCSA ratio
% y   :  Noisy  signal
% fs  : sampling frequency 
% Nh  : The desired Nh  number of negative eigenvaules 
% gm  : SCSA parameter
 
%% output parapeters
% h_denoising :  the optimal h 
% y_denoised   :  denoised  signal
% Nh_Dnz : Number of Eigen values needed fo denoising

% 
% function [h_full,Nh_full,y_full, Kh_op, K_Nh, y_denoised, h_denoising, Nh_Dnz ]=MRS_Denoising_KSCSA(K,gm,t,y,fs)%undesired_peak

function [h_full,Nh_full, y_full ]=Denoising_KSCSA(K,gm,t,y,fs)%undesired_peak
     global Plot_fig   

    Plot_fig=0;
    
    fprintf('\n  --> Full reconstruction of the noisy signal....')
    % [h, y_h, Nh,psinnor,kappa,Ymin,hlist0]=Full_reconstruct_SCSA_1D(y,fs,gm);
    [h_full, y_full, Nh_full]=Full_reconstruct_SCSA_1D_Dynamic_Scanning(y,fs,gm);
%     close all;figure; plot(y); hold on ;plot(y_full); hold on
    PSNR=psnr(y_full,y);
    fprintf('\n--> MRS denoising is fully reconstructed using  h=%f, PSNR=%d!!\n\n',h_full, PSNR)



%     Nh_dnz0=floor(Nh_full/3);
%     Plot_fig=1;
%     [Kh_op, h_denoising,  y_denoised,Nh_Dnz]= SCSA_1D_Nh_scanning(t,y,fs,Nh_dnz0,h_full,gm);
% 
% 
%     K_Nh=(Nh/Nh_Dnz);


function [K_op, h_op, yscsa,Nh_op]= SCSA_1D_Nh_scanning(t,y,fs,Nh0,h0, gm)
%     h=max(y);N=max(size(y));
global Plot_fig  shif figr

%% Apply h for full recostruction to get Nh

    [h0, yscsa,Nh_full]= SCSA_1D(y,fs,h0,gm);
    
    k=0;
    
    lgnd{1}=strcat('Noisy MRS signal');     lgnd{2}=strcat('Full reconstucted : h = ', num2str(h0) , ',   K_h = ', num2str(1), ',  Nh = ', num2str(Nh_full) , ',   K = ', num2str(1));
    
    if Plot_fig==1
        figure(figr);   plot(t,y,'LineWidth',2.5);hold on
                        plot(t,yscsa, 'LineWidth',2);hold on
                        legend(lgnd)

    end

    
    for K=1.25:0.25:1.75

        h=K*h0;  k=k+1;
        % Plot the results 
        [h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);
        Cost_function(k)=Nh;

        K_Nh=(Nh_full/Nh);

                    %% Plot the denoising results 
        if Plot_fig==1
            lgnd{k+2}=strcat('Denoised: h = ', num2str(h) , ',   K_h = ', num2str(K), ',  Nh = ', num2str(Nh) , ',   K = ', num2str(K_Nh));
            figure(figr);
                plot(t, yscsa + k*4 ,'LineWidth',2);hold on
                legend(lgnd);%({'Noisy input spectrum ', 'Denoised Spectrum ', 'Residual'},'Location','northwest');
                xlabel('t')
                ylabel('Intensity')
                set(gca,'YTickLabel',[])
                xlim([0 5])
                title(' SCSA MRS denoising K-SCSA');% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
                set(gcf,'color','w') 
                set(gca,'Xdir','reverse');
                set(gca,'fontsize',16)
            %     text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
                box 
            pause(0.3)

            
        end

         
         if Nh0>Nh
            K_op=K
            h_op=h;
            Nh_op=Nh;
            break;
         end
        
        
                 
    end

K_op=K   
h_op=h;
Nh_op=Nh;

 


function [h_full, y_full, Nh]=Full_reconstruct_SCSA_1D_Dynamic_Scanning(y,fs,gm)
global Plot_fig  
%     Plot_fig=1
    
if Plot_fig==1
    lgndr{1}='Noisy Signal';lgnd{1}='Noisy Signal';
    figure(1);
               subplot(211);plot(y);hold on ;legend(lgnd)
               subplot(212);plot(y);hold on ; legend(lgndr)
end

%% MRS signal 
y_positive=y-min(y);
h_max=0.3*(max(y_positive) + mean(y_positive));
% if max(y_positive)>1
    h_min=0.001;
% else
%     h_min=max(y_positive)/10;
% end

%% Setup the  Wait bar
global cnt_wait Tot_iter   wait_bar
Tot_iter=0;
nb_loops= 200;               % Number of loops 
M =2;                       % number of refinement loop to find the optimal h for denoising
for k=1:M
    Tot_iter= Tot_iter+floor(nb_loops/2^k);
end

cnt_wait=0;
wait_bar = waitbar(cnt_wait,'Searching for the optimal value of h .... ','Name','MRS signal denoising using SCSA');


%% Search for the optimal h 
eps=0.5;
for k=1:M
    h_list= linspace(h_min,h_max,floor(nb_loops/2^k));
    [h_op,Min_Cost]=search_4_optimal_h_reconstruction(y,fs,h_list,gm);
    h_min=h_op*(1-eps);
    h_max=h_op*(1+eps);
    %% plot the results     
    h_op_list(k)=h_op;
    Cost_function(k)=Min_Cost;
    
    
    % Plot the results 
    [h_op, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h_op,gm);
    Res_y=abs(yscsa-y);

    if Plot_fig==1
       lgnd{k+1}=strcat('h=',num2str(h_op),', Nh=',num2str(Nh));
       figure(100);  plot(Cost_function);hold off
%                      legend( 'MSE')
                     title(' Cost function maximization')

                     
        figure(1); subplot(211);plot((yscsa));title('Reconstructed signal') ;hold on ;legend(lgnd)

                   subplot(212);plot(Res_y);title(' Residual');hold on ; legend(lgnd)
                    %xlim([1 185])
        pause(0.3)
        d=1;
    end    
    
    
    
end


h_op=min(h_op_list(find(Cost_function==min(Cost_function))));

[h_full, y_full,Nh]= SCSA_1D(y,fs,h_op,gm);


d=1;

close(wait_bar)



function  [h_op,Min_Cost]=search_4_optimal_h_reconstruction(y,fs,h_list,gm)
global Tot_iter  wait_bar   cnt_wait
% start the loop for several values of h 
 for i=1:length(h_list)
 
    h = h_list(i);%0.045;%15.2;

    [h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA1D(y,fs,h,gm);
   
%     Th_y=0.98*max(y) ;
%     Maxy=abs(max(y)-max(yscsa));
%     if   Maxy <Th_y
%         Beta=0;
%     else
%         Beta=1e10;
%     end
    
    Cost_function(i)= mse(y,yscsa);%+ Beta*Maxy;

    

    fprintf('.')
    cnt_wait=cnt_wait+1;
    % Update waitbar and message
    waitbar(cnt_wait/Tot_iter,wait_bar)
%     figure(154);
%     plot(t,y,'LineWidth',2.5);hold on
%     plot(t, yscsa ,'LineWidth',1.2);hold off
% 
%     % hold on
%     % plot(f, y-yscsa-5,'g','LineWidth',1)
%     legend({'Noisy input spectrum ', 'Denoised Spectrum ','Residue'},'Location','northwest');
%     xlabel('t')
%     ylabel('Intensity')
%     set(gca,'YTickLabel',[])
%     xlim([0 5])
%     title([ ' SCSA MRS densoing with : h = ' num2str(h) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
%     set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
%     % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
%     box 
    
    
 end
Min_Cost=min(Cost_function);
h_op=min(h_list(find(Cost_function==min(Cost_function))));

% figure;plot(h_list,Cost_function)
d=1;





function [noise_areas,Peaks_areas]=find_metabolite_Noisy_ranges(y0,width_peaks)

y0=y0-min(y0);
%% Find tw peaks location 
N=max(size(y0));yf=y0;%(1:floor(0.48*N));
TH=mean(yf)+median(yf);
% yf_peaks=0*yf+TH;
% yf_peaks(find(yf>TH))=yf(find(yf>TH));
% Th=max(yf_smooth) /4;
% y_th=0*yf+TH;

idx=find(yf >TH);yf_peaks=0*yf;yf_peaks(idx)=yf(idx);

% close all; figure; plot(yf); hold on; plot(yf_peaks,'LineWidth',1.2); hold on;plot(y_th,'LineWidth',1.2); hold on; legend('yf','yf_peaks','y_th') 

yf_smooth   = sgolayfilt(yf_peaks, 2, 7);yf_smooth=yf_smooth-min(yf_smooth);

[pks,locs,w,p] = findpeaks(yf_peaks);

peaks_location=locs(find(yf_peaks(locs) > 0.7*mean(yf_peaks(locs))));
Peaks_areas=0;

Peak_TH=peaks_location(find( yf_smooth(peaks_location)>TH));

for i=1:max(size(Peak_TH))
    peaks_Locations(i,:)=[Peak_TH(i)-width_peaks,Peak_TH(i)+width_peaks];    
end

if find(peaks_Locations<=0)>0
    
    peaks_Locations=peaks_Locations(2:end,:);
end


for i=1:size(peaks_Locations,1)
   
    Peaks_areas=[Peaks_areas,peaks_Locations(i,1):peaks_Locations(i,2)];
    
end
Peaks_areas=unique(Peaks_areas(2:end));

indxes=1:max(size(y0));
noise_areas=setdiff(indxes,Peaks_areas);



function [yscsa, h_op, Nh]=SCSA_MRS_Denoising(SNR_regularization, t, yf, gm , fs , width_peaks)
% close all
global shif Plot_fig  

%% Get the noise and metabolites locations
[Noise,Metabolite]=find_metabolite_Noisy_Area(yf, width_peaks);

Noise(2)=Noise(1)+10;
%% Generate signals
%% MRS signal 
y = real(yf) ;%/max(real(yf))*76;

y_positive=y-min(y);
h_max=(max(y_positive) + mean(y_positive));
if max(y_positive)>1
    h_min=0.1;
else
    h_min=max(y_positive)/10;
end

fprintf('\n-->Searching for the optimal h. Please wait...')

% you can choose either a vecotr for h or one value 
% h_list=h_min : .1 :h_max ;  % search in this interval
% h_list=10:0.1:17;
eps=0.2;
% shif=0;close all;


%% Setup the  Wait bar
global cnt_wait Tot_iter   wait_bar
Tot_iter=0;
nb_loops= 90;               % Number of loops 
M =3;                       % number of refinement loop to find the optimal h for denoising
for k=1:M
    Tot_iter= Tot_iter+floor(nb_loops/2^k);
end

cnt_wait=0;
wait_bar = waitbar(cnt_wait,'Searching for the optimal value of h .... ','Name','MRS signal denoising using SCSA');


%% Search for the optimal h 
for k=1:M
    h_list= linspace(h_min,h_max,floor(nb_loops/2^k));
    [h_op,Min_Cost]=search_4_optimal_h_denoising(Metabolite, Noise,SNR_regularization, y,fs,h_list,gm);
    h_min=h_op*(1-eps);
    h_max=h_op*(1+eps);
    %% plot the results     
    h_op_list(k)=h_op;
    Cost_function(k)=Min_Cost;
    
end
h_op=min(h_op_list(find(Cost_function==min(Cost_function))));
[h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h_op,gm);

fprintf('\n--> MRS denoising is completed h=%f with Nh=%d!!\n\n',h_op,Nh )
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
        plot(t,y,'b','LineWidth',2.5);hold on
    end
    
    plot(t, yscsa+shif ,'r','LineWidth',2);hold on
    plot(t, y-yscsa-5,'g','LineWidth',1.5); hold off
    
    legend({'Noisy input spectrum ', 'Denoised Spectrum ', 'Residual'},'Location','northwest');
    xlabel('t')
    ylabel('Intensity')
    set(gca,'YTickLabel',[])
    xlim([0 5])
    title([ ' SCSA MRS denoising with : h = ' num2str(h_op) , '   Nh = ' num2str(Nh) , ',   \beta = ' num2str(SNR_regularization) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
    set(gcf,'color','w') 
    set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
    text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
    box 
end
    


function  [h_op,Min_Cost]=search_4_optimal_h_denoising(Metabolite, Noise,SNR_regularization, y,fs,h_list,gm)
global Tot_iter  wait_bar   cnt_wait
% start the loop for several values of h 
 for i=1:length(h_list)
 
    h = h_list(i);%0.045;%15.2;

    [Cost_SNR,Cost_peak,h, yscsa,Nh]=Cost_function_4_denoising(Metabolite, Noise,SNR_regularization, y,fs,h,gm);
    Cost_func_SNR(i)=Cost_SNR;
    Cost_func_MSE(i)=Cost_peak;
    Cost_function(i)=Cost_SNR+Cost_peak;

    fprintf('.')
    cnt_wait=cnt_wait+1;
    % Update waitbar and message
    waitbar(cnt_wait/Tot_iter,wait_bar)
%     figure(154);
%     plot(t,y,'LineWidth',2.5);hold on
%     plot(t, yscsa ,'LineWidth',1.2);hold off
% 
%     % hold on
%     % plot(f, y-yscsa-5,'g','LineWidth',1)
%     legend({'Noisy input spectrum ', 'Denoised Spectrum ','Residue'},'Location','northwest');
%     xlabel('t')
%     ylabel('Intensity')
%     set(gca,'YTickLabel',[])
%     xlim([0 5])
%     title([ ' SCSA MRS densoing with : h = ' num2str(h) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
%     set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
%     % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
%     box 
    
    
 end
Min_Cost=min(Cost_function);
h_op=min(h_list(find(Cost_function==min(Cost_function))));

function [noise_area,peaks_Locations]=find_metabolite_Noisy_Area(y0, width_peaks)
y0=y0-min(y0);
%% Find tw peaks location 
N=max(size(y0));yf=y0;%(1:floor(0.48*N));
TH=median(yf)+mean(yf(1:10));
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

Metabolite_area=1:peaks_location(end);

Peaks_areas=unique(Peaks_areas(2:end));
noise_areas=setdiff(Metabolite_area,Peaks_areas);
noise_area=[noise_areas(1),noise_areas(end)];


function [Cost_SNR,Cost_peak,h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]=Cost_function_4_denoising(Metabolite, Noise,SNR_regularization, y,fs,h,gm)

        % yscsa = scsa(h);
[h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);

% coco1(i) = norm(y(a1:b1) - yscsa(a1:b1));
% coco2(i) =  norm(y(a2:b2) - yscsa(a2:b2));
% coco3(i) =  norm(y(a3:b3) - yscsa(a3:b3));
% % coco(i) =norm(y - yscsa);

y_res=y-yscsa;

SNR_y=max(yscsa)/std(y(Noise(1):Noise(2)));   %second term retaed to the  SNR
SNR_scsa=max(yscsa)/std(yscsa(Noise(1):Noise(2)));   %second term retaed to the  residual SNR

Cost_SNR=(SNR_regularization/abs(SNR_scsa));%*abs(SNR_y);

Cost_peak=0;
for k=1:size(Metabolite,1) 
    Cost_peak=Cost_peak + norm(y_res(Metabolite(k,1:2)));
end


 