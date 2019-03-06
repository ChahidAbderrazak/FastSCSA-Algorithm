%% MRS signal Denosing  using SCSA + FFT: 
% This function denoises the real signal after projecting it in FFT

%% ######################  PARAMETERS ######################################
% Input 
% y       : Noisy real art of the  complex MRS FID signal
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


function [yf_scsa, h_op, fs_op, Nh]=SCSA_FFT_Denoising( yt, gm)

%% FFT Spectrum 
yf=fft(yt);

% close all;figure; plot(yt); hold on ; plot(abs(yf)); hold on; plot(angle(yf)); hold on

%% Get optimal denoising
[yf_scsa, h_op, fs_op, Nh]=SCSA_Denoising_RMSE( yf, gm);

%% IFFT Spectrum 
yf_scsa=ifft(yf_scsa);

d=1;