function [peak_pos, Big_peak_Loc]=find_peak_location(y)


y_smooth   = sgolayfilt(y, 1, 11); 
y_smooth   = sgolayfilt(y_smooth, 1, 21); 

% y_smooth=abs(zf_real_smooth); y_smooth=zf_real_smooth-min(y_smooth);

%% get the smoothed positive and negative parts of the signal
Mean_y=mean(y_smooth);
y_smooth_p=y_smooth*0;idx=find(y_smooth>=Mean_y);
y_smooth_p(idx)=y_smooth(idx);
y_smooth_n=y_smooth*0;idx=find(y_smooth<Mean_y);
y_smooth_n(idx)=-y_smooth(idx);

% close all; figure; plot(y_smooth_p); hold on; plot(y_smooth_n); hold off; 

%% Get the peaks location in from the positive and negative parts
peak_pos=[];
M_p=mean(y_smooth_p);Maxp= max(y_smooth_p);
if 1%abs(Maxp)>10*abs(M_p)
    peak_pos_p= Get_peak_pos(y_smooth_p);
    peak_pos=[peak_pos peak_pos_p];
end
% close all; figure; plot(y_smooth_p);

M_n=mean(y_smooth_n); Maxn=max(y_smooth_n) ;
if abs(Maxn)>5*abs(M_n)
    peak_pos_n= Get_peak_pos(y_smooth_n);
    peak_pos=[peak_pos peak_pos_n];

end

% close all; figure; plot(y_smooth_n);

d=1;

% peak_pos=[peak_pos_p peak_pos_n];
% close all; figure; plot(y); hold on ;plot(y_smooth); hold on ; vline(peak_pos)


peak_pos=sort(peak_pos);


a=0;
for k=1:max(size(peak_pos))-1
    peak=peak_pos(k)-3:peak_pos(k)+3;

    b=floor((peak_pos(k)+peak_pos(k+1))*0.5 );

    
     Min=min(y(peak_pos(k)));
     if Min>=0
         coef(k)=1;
     else
         coef(k)=-1;
     end
     area(k) = trapz(coef(k)*y(peak));
     
     peak_Loc(k,:)=[ a+1   b coef(k)]; a=b  ; 

end



% close all; figure; plot(y); hold on ;plot(y_smooth); hold on ; vline(peak_Loc(:,1))

  Min=min(y(peak_pos(k+1)));
     if Min>=0
         coef(k+1)=1;
     else
         coef(k+1)=-1;
     end
    peak=peak_pos(k)-3:peak_pos(k)+3;

 peak_Loc(k+1,:)=[ a+1   max(size(y)) coef(k+1)]; 
 peak=peak_pos(k+1)-3:peak_pos(k+1)+3;

 area(k+1) = trapz(coef(k+1)*y(peak));

 
 %% remove the small peaks
 idx=find(area>max(area)/5);

 
 
Big_peak_Loc=peak_Loc(idx,:);
Big_peak_Loc(1,1)=1;
Big_peak_Loc(end,2)=max(size(y));

peak_pos = Big_peak_Loc(:,1);
% peak_pos = peak_pos(idx);

% hold on; vline(peak_pos);
% vline(Big_peak_Loc(:,1));

d=1;



function peak_pos= Get_peak_pos(y)

N=max(size(y));
x=1:N;
[pks,locs,widths,proms] = findpeaks(y,x);
Peaks=sort(pks,'descend');
peaks_widths=sort(widths,'descend');
gain=Peaks.*widths;
% close all; figure; plot(y); hold on ;vline(locs); hold on ;   
% 
peaks_gain=sort(gain,'descend');
Max=max(y);
Np=max(size(Peaks));
cnt=1;
for k=1:2
    if Peaks(k)>Max/10
        Id0(cnt)=locs(find(pks==Peaks(k)));cnt=cnt+1;
        
    end
% Id1(k)=locs(find(widths==peaks_widths(k)))
% Id2(k)=locs(find(gain==peaks_gain(k)))

% peak_pos(k)=Id0

end
peak_pos=Id0;
% peak_pos=locs;

% POS=[Id0 Id1]
% [t,peak_pos]=element_frequency( POS)


function [t,peak_pos]=element_frequency(x)

xx = unique(x);       % temp vector of vals
x = sort(x);          % sorted input aligns with temp (lowest to highest)
t = zeros(size(xx)); % vector for freqs
% frequency for each value

k=1;
peak_pos=0;
for i = 1:length(xx)
    t(i) = sum(x == xx(i));
    
    if t(i)>1
        peak_pos(k)=xx(i);
        k=k+1;
    end
end

d=1;

