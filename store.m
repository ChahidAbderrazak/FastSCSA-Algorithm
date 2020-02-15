
    %% ##################  Denoising  MRS signal using  SCSA Cost function  ##################   
    fprintf('\n--> MSR Signal  Denoising SCSA using  Cost function  without refernce signal')
    shif=0;% close all
    SNR_regularization=500;  %  The regularization coefficient used for the optimization problem ~ 500;
    width_peaks=5;          %  The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
    figr=1;
    [yf_denoised,h_dnz,Nh_dnz]=SCSA_MRS_Denoising_regularization_no_ref(SNR_regularization, f, yf, gm , fs, width_peaks);
    % close all;figure; plot(yf); hold on ;plot(yf_denoised0); hold on;plot(yf_denoised); hold on

    %% ##################  Full reconstrucion of noisy   MRS signal  ##################  
    fprintf('\n--> Full reconstrucion of noisy   MRS signal ')
%     figr=2;
%     fprintf('\n--> MSR Signal  Denoising using K-SCSA')
    K=1.4;
    [h_full,Nh_full,y_full]=Denoising_KSCSA(K,gm,f,yf,fs);
    % close all;figure; plot(yf); hold on ; plot(yf_denoised); hold on; plot(y_full); hold on

     table_param(cnt,:)=[noise_level,h_full,h_dnz0,h_dnz,Nh_full,Nh_dnz0,Nh_dnz];

     cnt=cnt+1;