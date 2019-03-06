function [Peaks_area]=Get_Peak_Area_STD(y)
N=max(size(y));
Idx=find(y==max(y));
Peaks_area=Idx-3:Idx+3;
