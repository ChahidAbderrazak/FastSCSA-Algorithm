%% Water Suppression with interpolation and Scaling function All the signal 
% This script removes the water peak in the MRS signal based on selectoion
% of a subset of eigenfunction from the schrodinger operaotor using an
% optimal paramter h that id defined using a adaptive interval scanning
%
% Important: the input data are stored in "./Input_data/*.mat' files  containning:
% [img : reference WS-MRS signal] , [original_file : name of source file] 
% [noisy_img :  MRS signal], [noisy_file : name of source file] 

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: Sep,  2017
%  
%% ###########################################################################
close all ;  clear all; clc; warning ('off','all');; tic
global  Results_path post_save_tag name_data  store_decomposition ID x_i_list scaling_EN_list  frame_EN_list TypeScal_list pulse_nb_list factor_list 
addpath ./Solver_Function ;Include_function ;ID=120;log_html_file; 
global Plot_fig  shif  figr;
Plot_fig=0; figr=1;

%% Set the results directory
Root='./Results/XXXX/';
Results_path=strcat(Root,'/',char(datetime('today')));mkdir(Results_path) % The obtained results figures  will be saved in:

%% load MRS spectrum
% Load_prepare_MRS_Spectum2018  %Load Inivo data sent on 18-03-31
cnt=1;

%% ################ Code start here  ######################
% SCSA Parameter
gm=0.5;fs=1;

for noise_level=10:5:30
    
    [f, y0_complex,y_complex ]=generate_Lorntzian(noise_level);
    yf=real(y_complex);    yf0=real(y0_complex);
    %% ##################  Denoising  MRS signal using  SCSA Cost function  ##################  
    fprintf('\n--> MSR Signal  Denoising SCSA using  PSNR ')
    [yf_denoised0,h_dnz0,Nh_dnz0]=SCSA_MRS_Denoising_PSNR(  f, yf0, yf, gm , fs );
    % close all;figure; plot(yf); hold on ; plot(yf_denoised0); hold on

    %% ##################  Denoising  MRS signal using  SCSA Cost function  ##################   
    fprintf('\n--> MSR Signal  Denoising SCSA using  Cost function ')
    shif=0;% close all
    SNR_regularization=500;  %  The regularization coefficient used for the optimization problem ~ 500;
    width_peaks=5;          %  The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
    figr=1;
    [yf_denoised,h_dnz,Nh_dnz]=SCSA_MRS_Denoising(SNR_regularization, f, yf, gm , fs, width_peaks);
    % close all;figure; plot(yf); hold on ;plot(yf_denoised0); hold on;plot(yf_denoised); hold on

    %% ##################  Full reconstrucion of noisy   MRS signal  ##################  
    fprintf('\n--> Full reconstrucion of noisy   MRS signal ')

%     figr=2;
%     fprintf('\n--> MSR Signal  Denoising using K-SCSA')
    K=1.4;
    [h_full,Nh_full, y_full,Kh_op, K_Nh, yf_denoised_K,Nh_dnzK]=MRS_Denoising_KSCSA(K,gm,f,yf,fs);
     % close all;figure; plot(yf); hold on ;plot(yf_denoised); hold on ;plot(yf_denoised_K); hold on

    [h_full,Nh_full,y_full]=Denoising_KSCSA(K,gm,f,yf,fs);
    % close all;figure; plot(yf); hold on ; plot(yf_denoised); hold on; plot(y_full); hold on

     table_param(cnt,:)=[noise_level,h_full,h_dnz0,h_dnz,Nh_full,Nh_dnz0,Nh_dnz];

     cnt=cnt+1;
end

%% plot results 
plot_NhRatio_Comparison(f,yf0, yf, y_full, yf_denoised, yf_denoised0, noise_level,h_full,h_dnz0,h_dnz,Nh_full,Nh_dnz0,Nh_dnz)    
    
%% Save  results 
Nh_Ratio=table_param(:,5)./table_param(:,6);
table_param=[table_param Nh_Ratio];
%% Add the last results
colnames={'Noise_Level','h_reconstruction','h_PSNR','h_previous','Nh_reconstruction','Nh_PSNR','Nh_previous', 'Nh_Ratio'};
perform_output= array2table(table_param, 'VariableNames',colnames)
save(strcat(Results_path,'/Nh_ratio_with_Noise.mat'),'perform_output','yf','yf0','h_full','h_dnz0','h_dnz')

%  %% Save the results 
%  Save_optimal_solution_Denoising_final_reduced  
%  
% %% THE END
% fprintf('\n\n THE END time=%0.2f  Sec\n',toc) ;
% 

d=1;