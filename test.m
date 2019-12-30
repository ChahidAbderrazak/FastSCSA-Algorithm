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

for noise_level=20%0:5:30
    
    %% Generate noisy academic or loratzian signals
    Generate_noisy_signal
%     load('.\input\pulse_signal.mat')

    %% ##################  Denoising  MRS signal using  SCSA Cost function  ##################  

    fprintf('\n-->  Signal  Denoising SCSA  using fs scanning. Cost function= PSNR')
    [yscsa0, h_op0, fs_op0, Nh_op0]=SCSA_Denoising_PSNR( y, gm , fs ,y0);
    PSNR_op0=psnr(y0,yscsa0)
    param0= strcat(' : PSNR =>[ h^*=',num2str(h_op0),' , f_s^*=',num2str(fs_op0),', Nh= ',num2str(Nh_op0) ,', PSNR=',num2str(PSNR_op0),']');
    
%     
%     fprintf('\n-->  Signal  Denoising SCSA  using fs scanning. Cost function= Smoothness')
%     [yscsa1, h_op1, fs_op1, Nh_op1,Noise_area1]=SCSA_Denoising_Smoothness( y, gm , fs);
%     PSNR_op1=psnr(y0,yscsa1)
%     param1= strcat(': Smoothness => [ h^*=',num2str(h_op1),' , f_s^*=',num2str(fs_op1),', Nh= ',num2str(Nh_op1) ,', PSNR=',num2str(PSNR_op1),']');
% 
%     
%     fprintf('\n-->  Signal  Denoising SCSA  using h . Cost function=  soft-thresolding [ smoothness + PRR]')
%     [yscsa2, h_op2, fs_op2, Nh_op2,Noise_area2]=SCSA_Denoising_PRR_smoothness(y, y, gm , fs);
%     PSNR_op2=psnr(y0,yscsa2)
%     param2= strcat(': Smoothness + PRR => [ h^*=',num2str(h_op2),' , f_s^*=',num2str(fs_op2),', Nh= ',num2str(Nh_op2) ,', PSNR=',num2str(PSNR_op2),']');
% 
%     
%     fprintf('\n-->  Signal  Denoising SCSA  using fs scanning. Cost function= Total area')
% 
%     [yscsa3, h_op3, fs_op3, Nh_op3]=SCSA_Denoising_Total_Area( y, gm , fs);
%     PSNR_op3=psnr(y0,yscsa3)
%     param3= strcat(': Total area => [ h^*=',num2str(h_op3),' , f_s^*=',num2str(fs_op3),', Nh= ',num2str(Nh_op3) ,', PSNR=',num2str(PSNR_op3),']');



    fprintf('\n-->  Signal  Denoising SCSA  using fs scanning. Cost function= RMSE')
    [yscsa4, h_op4, fs_op4, Nh_op4]=SCSA_Denoising_RMSE( y, gm);
    PSNR0=psnr(y0,y)
    PSNR_op4=psnr(y0,yscsa4)

    param4= strcat(' : RMSE =>[ h^*=',num2str(h_op4),' , f_s^*=',num2str(fs_op4),', Nh= ',num2str(Nh_op4) ,', PSNR=',num2str(PSNR_op4),']');



    
    %% plot the comprision
    Plot_denoisiong_Comparison(noise_level, t, y0, y,  yscsa0, param0, yscsa4, param4);
    
%     Plot_denoisiong_Comparison(noise_level, f_ppm, y0, yf,  yscsa0, param0, yscsa3, param3);
%     
%     
%     PSNR_op0
%     PSNR_op1
%     PSNR_op2
%     PSNR_op3
%     PSNR_op4
end


% %% Save for comparison with WATV
% 
% save('G:\My Drive\MyCodes\GitHub\InProcess_Projects\Methods_comparision\Signal denoising\2015 WATV\code\input\SCSA_results0.mat','yscsa0')