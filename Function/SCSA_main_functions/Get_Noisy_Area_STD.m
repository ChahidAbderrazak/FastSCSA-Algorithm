function [noise_areas,Peaks_area]=Get_Noisy_Area_STD(y0,STD_TH_Coef)

N=max(size(y0));yf=y0;%(1:floor(0.48*N));
W=floor(N/50);  %the sliding windows size
for k=1:N-W
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
