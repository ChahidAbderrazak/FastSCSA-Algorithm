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
close all ;  clear all; clc; warning ('off','all'); addpath ./Function ;Include_function; 
global Plot_fig  shif  figr h0 Results_path;

%% Input parameters
% SCSA Parameter
gm=0.5;fs=1;

% Set the results directory
Root='./Results/Eusipco2020/'; Results_path=strcat(Root,'/',char(datetime('today')));mkdir(Results_path) % The obtained results figures  will be saved in:
Plot_fig=0; figr=1;

%% ################ Code start here  ######################
fprintf('\n------------------------------------------ ')
fprintf('\n|  Signal  reconstruction using fastSCSA  | ')
fprintf('\n------------------------------------------ ')

%% load MRS spectrum
% Load_prepare_MRS_Spectum2018  %Load Inivo data sent on 18-03-31
cnt=1;
Noise_list='_Nz';

for sig_num= [-3:-1 3 6 7 ]
    
    for noise_level=10%:10:30

        Data_preparation_Load
        y_desired=yf;    type_estimation='Recons';    % Reconstruction 
        y_desired=yf0; type_estimation='Denoiz';    % Denoising 
        tag_experimt=strcat(type_estimation,'_',sig_name);
        Noise_list=strcat(Noise_list,'-',num2str(noise_level));
        used_signl{cnt}=sig_name(1:end-1);
        %% search interval 
        h_min=sum(sqrt(yf-min(yf)))/(pi*N);   %0.001;
        h_max=sum(sqrt(yf-min(yf)))/(pi);     %(max(yf) + mean(yf));

        %% Fast search
        fprintf('\n-->  fastSCSA ')

        [yf_fast,h_fast,Nh_fast, Nb_iter_fast]=FastSCSA_reconstruction_PSNR(f, y_desired, yf, gm , fs ,h_min, h_max);
        MSE_fast=mse(yf_fast,y_desired)

        %% grid search
        fprintf('\n-->  Grid search ')

        h_step=0.1;
        h0=h_fast ; % just to accelerate the comparison
        [yf_grid, h_grid, Nh_grid, Nb_iter_grid]=SCSA_reconstruction_scanning(y_desired, yf, gm , fs ,h_min, h_max, h_step);
        MSE_grid=mse(yf_grid,y_desired)

        %% Algorithm Comparison
         table_param(cnt,:)= [noise_level, h_grid ,h_fast, Nh_grid, Nh_fast, MSE_grid, MSE_fast, Nb_iter_grid ,Nb_iter_fast ];
         cnt=cnt+1;



        % plot figure
        plot_all_reconstructed_signals

    end

end
d=1;


%% Save  results 
colnames={'Noise_Level','h_grid' ,'h_fast','Nh_grid','Nh_fast','MSE_grid','MSE_fast','iter_grid','iter_fast'};
perform_output= array2table(table_param, 'VariableNames',colnames);
perform_output.Singal=used_signl'
name_report=strcat('/',tag_experimt,'_report',Noise_list,'.csv');
writetable(perform_output,strcat(Results_path,name_report))

