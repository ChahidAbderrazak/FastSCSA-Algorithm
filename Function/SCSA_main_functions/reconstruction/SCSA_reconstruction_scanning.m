%% MRS signal Denosing  using SCSA Method based on PSNR: 
% This function search for  optimal value of the h the ensures
% better estimatin of y_in from y_desired by minimazing the error  
%  |yscsa-y_in|
%% ######################  PARAMETERS ######################################
% Input 
% t         : The MRS signal fequancies in t
% y_desired : desired estimation 
% y      : input  signal
% gm, fs  : SCSA parameter: gm=0.5, fs=1
% h_min, h_max: initial search interval 

% Output
% yscsa: The reconstructed  signal
% h_op : The optimal paramter h
% Nh_op   : The number of negative eigenvalues used for the reconstruction

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: June,  2018
% King Abdullah University of Sciences and Technology (KAUST)


function [yscsa, h_op, Nh_op,Nb_iteration]=SCSA_reconstruction_scanning(y_desired, y_in, gm , fs ,h_min, h_max, h_step)
global h0
L=200;
mse_op=inf;
h_list= h_min:h_step:h_max;
N=length(h_list);
idx=find(h_list<=h0);
idx0=idx(end);

if idx0+L>N; Idmax=N; else;Idmax=idx0+L;end
if idx0<L;Idmin=1;else;Idmin=idx0-L;end

h_list=h_list(Idmin:Idmax);

for i=1:length(h_list)
    tic
    h = h_list(i);%0.045;%15.2;
    tic
    [h, y_n,Nh]= SCSA_1D(y_in,fs,h,gm);
    tim_unit=toc;
    mse_i=mse(y_desired,y_n);
    cost_msei(i)=mse_i;
    if mse_op>mse_i
        h_op=h
        Nh_op=Nh;
        mse_op=mse_i;
        yscsa=y_n;
    end
end

Nb_iteration=length(h_list);
figure;
plot(h_list,cost_msei)