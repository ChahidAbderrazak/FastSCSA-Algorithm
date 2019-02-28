%%  SCSA pulse shaped signal denosing 
% This script  denoise real values pulse shaped signals

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
% Done: Mar,  2019
%  
%% ###########################################################################
close all ;  clear all; clc; tic
global  Results_path post_save_tag name_data  store_decomposition ID x_i_list scaling_EN_list  frame_EN_list TypeScal_list pulse_nb_list factor_list 
addpath ./Function ;Include_function ;

%% Set the results directory
Root='r:/SCSA_Algorithm_Matla/Understand_SCSA/';
Results_path=strcat(Root,'Results/',char(datetime('today'))); % The obtained results figures  will be saved in:

%% ################ Code start here  ######################
gm=0.5;fs=1;

for noise_level=10%0:5:30
    
    [f_ppm,zt0,zt,zf0,zf ]=generate_Lorntzian(noise_level);
    yf=real(zf);    yf0=real(zf0);
    % close all;figure; plot(yf); hold on ; plot(yf0); hold on
    yt=real(zt);    yt0=real(zt0);
    % close all;figure; plot(yt); hold on ; plot(yt0); hold on

    %% ##################  Denoising  MRS signal using  SCSA Cost function  ##################  

    fprintf('\n-->  Signal  Denoising SCSA  using fs scanning')
    [yscsa1, h_op1, fs_op1, Nh,Noise_area]=SCSA_Denoising_Scanning( yf, gm , fs ,yf0);
    PNSR_op1=psnr(yf,yscsa1);

    
    fprintf('\n-->  Signal  Denoising SCSA  using h soft-thresolding')

    [yscsa2, h_op2, Nh,Noise_area]=SCSA_Denoising_PRR_smoothness(yf, yf, gm , fs);
    PNSR_op2=psnr(yf,yf_denoised');

end