clc;
clear 

SNR_dB =     0  ;
SNR = SNR_dB;
SNR_linear=10.^(SNR/10.);
sigma2=1/SNR_linear;

nbrOfSetups =  1000;

%%% system parameters
N =   100:40:300;%   128:64:384; % number of beams (transmit antennas)
len = length(N);

K = 1;
L = 10; % number of all paths
kappa = 10;
sparsity_f = 8 ;
sparsity_f_1 = 4 ;

u = 0.8;
sparsity_1 = 1;

% us = 0.04;
% sparsity_1s = 0.6;
sparsity_1s = [0.4,       0.4,       0.4,       0.4,      0.4,       0.4  ];
us =               [ 0.08,     0.08,     0.08,     0.08,    0.08,     0.08] ;

gamma  =  5/10  ;
Lf =floor( L*gamma); % number of paths for far-field
Ln = floor(L*(1-gamma)); % number of paths for near-field

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space


NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );
NMSE_FF_OMP = zeros(1,len);
NMSE_FF_SGP = zeros(1,len );
NMSE_NF_OMP = zeros(1,len);
NMSE_NF_SGP = zeros(1,len );
NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

NMSE_HF_SGP_2 = zeros(1,len );
NMSE_HF_SGP_3 = zeros(1,len );
% NMSE_HF_SGP_4 = zeros(1,len );
% NMSE_HF_SGP_5 = zeros(1,len );
% NMSE_HF_SGP_6 = zeros(1,len );

z=0;
t0 = clock;

for  i = 1 : len

    DFT  = (1/sqrt(N(i)))*exp(-1i*pi*[0:N(i)-1]'*[-(N(i)-1)/2:1:(N(i)/2)]*(2/N(i))) ;
    label_theta = [-(N(i)-1)/2:1:(N(i)/2)]*(2/N(i));

    % the near-field polar-domain transform matrix [5]
    Rmin=10;
    Rmax=80;
    eta = 2.5; 
    [Polarcodebook, label, dict_cell, label_cell] = QuaCode(N(i), d, lambda_c, eta, Rmin, Rmax);
    % [Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
    S = size(Polarcodebook, 2);

    H = zeros(N(i),K) ;
    % Y = zeros(M,K) ;
    Y = zeros(N(i),K) ;

%     Rh=zeros(N(i) , N(i) );
%     for s=1:5000
%         [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N(i), Lf, Ln, d, fc,Rmin, Rmax, kappa);
%         Rh=Rh+h*h';
%     end
%     Rh=Rh./(5000);

    error_HF_OMP = 0;  
    error_HF_SGP = 0;     
%     error_FF_OMP = 0;  
%     error_FF_SGP = 0;     
%     error_NF_OMP = 0;  
%     error_NF_SGP = 0;             
    error_LS = 0;  
    error_MMSE = 0;
    error_HF_OMP_1 = 0;  
    error_HF_SGP_1 = 0; 
    error_HF_SGP_2 = 0; 
    error_HF_SGP_3 = 0; 
%     error_HF_SGP_4 = 0; 
%     error_HF_SGP_5 = 0; 
%     error_HF_SGP_6 = 0; 

    error_FF_OMP_1 = 0;

    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        for k = 1 : K
            [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N(i), Lf, Ln, d, fc,Rmin, Rmax, kappa);
            noise = sqrt(sigma2)*(randn(N(i),1)+1i*randn(N(i),1))/sqrt(2);

            y=  h+ noise ;
            H(:,k) = h;
            Y(:,k) = y;
        end

z = z +  norm(h)/norm(noise)/     nbrOfSetups;



        %% the hybrid-field OMP based scheme
        Hhat_homp = Hybrid_OMP( Y , DFT, Polarcodebook,    Lf*sparsity_f ,   Ln*sparsity_f  );
        error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;
        
        %% the hybrid-field OMP 1 1 based scheme
        Hhat_homp_1 = Hybrid_OMP_1( Y , DFT     , Polarcodebook    ,    L *sparsity_f_1  ,L);
        error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;


        %% the proposed hybrid-field SGP based scheme
        
        [  Hhat_HF_SGP  , support , g ,sup_f ,sup_n   ]=Hybrid_SGP( Y, DFT, Polarcodebook,   Lf  *sparsity_1 ,    Ln *sparsity_1 , u );

        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

%         [A, G       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 0.1, 0.1, 1e-8 );
%         error_HF_SGP_3 =  error_HF_SGP_3 +   norm(    H -  A*G   ,'fro')^2/norm( H   ,'fro')^2  ;

        %% the proposed hybrid-field SGP based scheme
        
        [  Hhat_HF_SGP_1  , support , g ,sup_f ,sup_n  ]=Hybrid_SGP_1( Y, DFT, Polarcodebook,   L *sparsity_1s(i) , us(i) ,   L);

        error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;

        [A, G       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 0.05, 0.005, 1e-8 );
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H -  A*G   ,'fro')^2/norm( H   ,'fro')^2  ;

%         [A, G       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 0.1, 0.01, 1e-8 );
%         error_HF_SGP_3 =  error_HF_SGP_3 +   norm(    H -  A*G   ,'fro')^2/norm( H   ,'fro')^2  ;        


%         %% the LS
%         Hhat_LS   =  h + sqrt(sigma2)*(randn(N(i),1)+1i*randn(N(i),1))/sqrt(2);
% %         Hhat_LS   =   Y  ;
% 
%         error_LS    =  error_LS +     norm(    H -  Hhat_LS  ,  'fro'   )^2   /    norm( H  ,'fro')^2  ;

%        %% the MMSE
%        hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2*eye(N(i))  )  )  ) * Hhat_LS; 
% 
%        error_MMSE =  error_MMSE    +     norm(    H -  hhat_MMSE  ,  'fro')^2/norm( H  ,'fro')^2  ;      

    end




NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
NMSE_LS(i) = error_LS / nbrOfSetups ;
NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;
NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;
NMSE_HF_SGP_2(i) = error_HF_SGP_2 /  nbrOfSetups ;
NMSE_HF_SGP_3(i) = error_HF_SGP_3 /  nbrOfSetups ;

end





NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;
NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;
NMSE_HF_SGP_2 = 10*log10(NMSE_HF_SGP_2)   ;
NMSE_HF_SGP_3 = 10*log10(NMSE_HF_SGP_3)   ;

figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  N,   NMSE_LS               ,      '     k   -  o  '    , 'linewidth',  1);
% plot(  N,   NMSE_MMSE               ,      '     k   -  ^  '    , 'linewidth',  1);
plot(  N,   NMSE_HF_OMP               ,      '    b   -  s  '    , 'linewidth',  1.5);
plot(  N,   NMSE_HF_OMP_1               ,      '     b   --  o  '    , 'linewidth',  1.5);
plot(  N,   NMSE_HF_SGP   ,      '     r  -   p '    ,  'linewidth' , 1.5);
% plot(  N,   NMSE_HF_SGP_1   ,      '     r  --   d '    ,  'linewidth' , 1);
plot(  N,   NMSE_HF_SGP_2   ,      '     r  --   d '    ,  'linewidth' , 1.5);
% plot(  N,   NMSE_HF_SGP_3   ,      '     g  --   o '    ,  'linewidth' , 2);
xlabel('Number of antennas at the BS','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
legend( 'Hybrid-field OMP (with $\gamma$) [23]'  ,'Hybrid-field OMP (without $\gamma$) [27]'  ,  'On-grid hybrid-field SGP (with $\gamma$)' ,  'Off-grid hybrid-field SGP (without $\gamma$)', 'Interpreter','Latex');
% xlim([128 512])
ylim([-8 1])

% figure('color',[1,1,1]); hold on; box on; grid on;
% % plot(  N,   NMSE_LS               ,      '     k   -  o  '    , 'linewidth',  1);
% plot(  N(3:1:7),   NMSE_MMSE(3:1:7)               ,      '     b   -  ^  '    , 'linewidth',  1);
% plot(  N(3:1:7),   NMSE_HF_OMP(3:1:7)               ,      '     r   -  s  '    , 'linewidth',  1);
% plot(  N(3:1:7),   NMSE_HF_OMP_1(3:1:7)               ,      '     r   --  s  '    , 'linewidth',  1);
% plot(  N(3:1:7),   NMSE_HF_SGP(3:1:7)   ,      '     g  -   p '    ,  'linewidth' , 1);
% plot(  N(3:1:7),   NMSE_HF_SGP_1(3:1:7)   ,      '     g  --   p '    ,  'linewidth' , 1);
% xlabel('The number of antennas at the BS','Interpreter','Latex');
% ylabel('NMSE (dB)','Interpreter','Latex');
% legend(  'MMSE','hybrid-field-OMP',    'Improved hybrid-field-OMP','hybrid-field-SGP'  ,     'Improved hybrid-field-SGP', 'Interpreter','Latex');
% xlim([256 512])
% ylim([-4 4.5])


