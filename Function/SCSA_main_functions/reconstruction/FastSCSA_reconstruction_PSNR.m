%% MRS signal Denosing  using SCSA Method based on PSNR: 
% This function search for  optimal value of the h the ensures
% better estimatin of y_in from y_desired by minimazing the error  
%  |yscsa-y_in|
%% ######################  PARAMETERS ######################################
% Input 
% t         : The MRS signal fequancies in t
% y_desired : desired estimation 
% y_in      : input  signal
% gm, fs  : SCSA parameter: gm=0.5, fs=1
% h_min, h_max: initial search interval 

% Output
% yscsa: The reconstructed  signal
% h_op : The optimal paramter h
% Nh   : The number of negative eigenvalues used for the reconstruction

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: June,  2018
% King Abdullah University of Sciences and Technology (KAUST)

function [yscsa, h_op, Nh, Nb_itera]=FastSCSA_reconstruction_PSNR(t, y_desired, y_in, gm , fs ,h_min, h_max)
% close all
global shif Plot_fig  Nb_itera 
%  input parameters
nb_loops= 90;               % Number of loops 
M =3;                       % number of refinement loop to find the optimal h for denoising
eps=0.2;                    % size of interval to search in around an optimal h 
% h_min=0.001;

%% sigal preparation
y = real(y_in) ;%/max(real(y_in))*76;
y0 = real(y_desired) ;%/max(real(y_in))*76;
% y_positive=y-min(y);
% h_max=(max(y_positive) + mean(y_positive));

fprintf('\n-->Searching for the optimal h using PSNR. Please wait...')

% you can choose either a vecotr for h or one value 
% h_list=h_min : .1 :h_max ;  % search in this interval
% h_list=10:0.1:17;
% shif=0;close all;


%% Setup the  Wait bar
global cnt_wait Tot_iter   wait_bar
Tot_iter=0;Nb_itera=0;

for k=1:M
    Tot_iter= Tot_iter+floor(nb_loops/2^k);
end

cnt_wait=0;
wait_bar = waitbar(cnt_wait,'Searching for the optimal value of h using PSNR error.... ','Name','Signal reconstruction FastSCSA');


%% Search for the optimal h 
for k=1:M
    h_list= linspace(h_min,h_max,floor(nb_loops/2^k));
    [h_op,Min_Cost]=search_4_optimal_h_denoising_PSNR( y0, y,fs,h_list,gm);
    h_min=h_op*(1-eps);
    h_max=h_op*(1+eps);
    %% plot the results     
    h_op_list(k)=h_op;
    Cost_function(k)=Min_Cost;
    
end
h_op=min(h_op_list(find(Cost_function==min(Cost_function))));
[h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h_op,gm);

PSNR=psnr(yscsa,y0);
fprintf('\n--> Signal reconstruction is completed h=%f, PSNR=%d!!\n\n',h_op, PSNR)
close(wait_bar)


 %% Get the optimal     
% figure(2)
% plot(Cost_function); hold on
% plot(Cost_func_SNR); hold on
% plot(Cost_func_MSE); hold on
% 
% xlabel('Iterations')
% ylabel('Cost function ')
% title([ 'The Iterative MRS signal denoising']);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])


%% Plot the denoising results 
if Plot_fig==1
      figure;
    if shif==0
        plot(t,y0,'k','LineWidth',2.5);hold on
        plot(t,y,'b','LineWidth',2.5);hold on
    end
    
    plot(t, yscsa ,'r','LineWidth',2);hold on
    plot(t, y-yscsa-(0.45*min(yscsa)),'g','LineWidth',1.5); hold off
    
    legend({'reference input signal ', 'input signal ', 'Estimated signal ', 'Residual'},'Location','northeast');
    xlabel('t')
    ylabel('Intensity')
    set(gca,'YTickLabel',[],'XTickLabel',[])
%     xlim([t(1) t(end)])
    title([ ' SCSA signal estimation with : h = ' num2str(h_op) , '   Nh = ' num2str(Nh) , ',   PSNR = ' num2str(PSNR) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
    set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
%     text(1.8,60,'NAA');text(1.40,40,'Lac(1)');text(1.2,30,'Lac(2)');text(4.7,30,'Water residue');
    box 
end
    


function  [h_op,Min_Cost]=search_4_optimal_h_denoising_PSNR(y0, y,fs,h_list,gm)
global Tot_iter  wait_bar   cnt_wait Nb_itera
% start the loop for several values of h 
 for i=1:length(h_list)
 
    h = h_list(i);%0.045;%15.2;
    [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);
    Nb_itera=Nb_itera+1;
    Cost_function(i)= psnr(y0,yscsa);

    fprintf('.')
    cnt_wait=cnt_wait+1;
    % Update waitbar and message
    waitbar(cnt_wait/Tot_iter,wait_bar)
%     figure(154);
%     plot(t,y,'LineWidth',2.5);hold on
%     plot(t, yscsa ,'LineWidth',1.2);hold off
% 
%     % hold on
%     % plot(f, y-yscsa-5,'g','LineWidth',1)
%     legend({'Noisy input signal ', 'Denoised signal ','Residue'},'Location','northwest');
%     xlabel('t')
%     ylabel('Intensity')
%     set(gca,'YTickLabel',[])
%     xlim([0 5])
%     title([ ' SCSA MRS densoing with : h = ' num2str(h) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
%     set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
%     % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
%     box 
    
    
 end
Min_Cost=max(Cost_function);
h_op=min(h_list(find(Cost_function==max(Cost_function))));

% close all; figure; plot(Cost_function)
d=1;



function [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm)

Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));
N=max(size(y));
%% remove the negative part
Ymin=min(y);

 y_scsa = y -Ymin;
%% Build Delta metrix for the SC_hSA
feh = 2*pi/N;
D=delta(N,fs,feh);

%% start the SC_hSA
Y = diag(y_scsa);
SC_h = -h*h*D-Y; % The Schrodinger operaor

% = = = = = = Begin : The eigenvalues and eigenfunctions
[psi,lamda] = eig(SC_h); % All eigenvalues and associated eigenfunction of the schrodinger operator
% Choosing  eigenvalues
All_lamda = diag(lamda);
ind = find(All_lamda<0);


%  negative eigenvalues
Neg_lamda = All_lamda(ind);
kappa = diag((abs(Neg_lamda)).^gm); 
Nh = size(kappa,1); %%#ok<NASGU> % number of negative eigenvalues



if Nh~=0
    
% Associated eigenfunction and normalization
psin = psi(:,ind(:,1)); % The associated eigenfunction of the negarive eigenvalues
I = simp(psin.^2,fs); % Normalization of the eigenfunction 
psinnor = psin./sqrt(I);  % The L^2 normalized eigenfunction 


%yscsa =4*h*sum((psinnor.^2)*kappa,2); % The 1D SC_hSA formula
yscsa1 =((h/Lcl)*sum((psinnor.^2)*kappa,2)).^(2/(1+2*gm));
else
    
  psinnor = 0*psi;  % The L^2 normalized eigenfunction 

  yscsa1=0*y;
  yscsa1=yscsa1-10*abs(max(y));
  disp('There are no negative eigenvalues. Please change the SCSA parameters: h, gm ')
end


if size(y_scsa) ~= size(yscsa1)
yscsa1 = yscsa1';
end
 
 %% add the removed negative part
 yscsa = yscsa1 + Ymin;


squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;




    %**********************************************************************
    %*********              Numerical integration                 *********
    %**********************************************************************

    % Author: Taous Meriem Laleg

    function y = simp(f,dx);
    %  This function returns the numerical integration of a function f
    %  using the Simpson method

    n=length(f);
    I(1)=1/3*f(1)*dx;
    I(2)=1/3*(f(1)+f(2))*dx;

    for i=3:n
        if(mod(i,2)==0)
            I(i)=I(i-1)+(1/3*f(i)+1/3*f(i-1))*dx;
        else
            I(i)=I(i-1)+(1/3*f(i)+f(i-1))*dx;
        end
    end
    y=I(n);
    

    %**********************************************************************
    %*********             Delata Metrix discretization           *********
    %**********************************************************************
    
    
%Author: Zineb Kaisserli

function [Dx]=delta(n,fex,feh)
    ex = kron([(n-1):-1:1],ones(n,1));
    if mod(n,2)==0
        dx = -pi^2/(3*feh^2)-(1/6)*ones(n,1);
        test_bx = -(-1).^ex*(0.5)./(sin(ex*feh*0.5).^2);
        test_tx =  -(-1).^(-ex)*(0.5)./(sin((-ex)*feh*0.5).^2);
    else
        dx = -pi^2/(3*feh^2)-(1/12)*ones(n,1);
        test_bx = -0.5*((-1).^ex).*cot(ex*feh*0.5)./(sin(ex*feh*0.5));
        test_tx = -0.5*((-1).^(-ex)).*cot((-ex)*feh*0.5)./(sin((-ex)*feh*0.5));
    end
    Ex = full(spdiags([test_bx dx test_tx],[-(n-1):0 (n-1):-1:1],n,n));
    
    Dx=(feh/fex)^2*Ex;

    



