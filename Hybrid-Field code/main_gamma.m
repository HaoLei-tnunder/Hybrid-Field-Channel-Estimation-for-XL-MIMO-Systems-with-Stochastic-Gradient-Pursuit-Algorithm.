clc;
clear 


SNR_dB =    0 ;

SNR = SNR_dB;
SNR_linear=10.^(SNR/10.);
sigma2=1/SNR_linear;

nbrOfSetups = 1000;

%%% system parameters
N = 256; % number of beams (transmit antennas)
K = 1;
L = 10; % number of all paths
kappa = 10;

sparsity_f = 8 ;
sparsity_f_1 = 4 ;

u = 0.8;
sparsity= 1;
us = 0.08;
sparsity_s= 0.4;


gamma  =  0: (2/L)  :1 ;

len =   length(  gamma  );  

% M = 256; % number of pilot overhead

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space

% the far-field angle-domain DFT matrix
% D = N; %字典规模
% row = (-(N - 1)/2:(N - 1)/2)' ;
% col = -1 + 2/D : 2/D : 1 ;
% DFT  =  exp( 1j*  pi * row * col ) / sqrt(N);

DFT  = (1/sqrt(N))*exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N)) ;
label_theta = [-(N-1)/2:1:(N/2)]*(2/N);

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5; 
[Polarcodebook, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
% [Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
S = size(Polarcodebook, 2);


NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );
NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_SGP_2 = zeros(1,len );
NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

H = zeros(N,K) ;
Y = zeros(N,K) ;

t0 = clock;

for  i = 1 :  len

    Lf = L*gamma(i); % number of paths for far-field
    Ln =floor( L*(1-gamma(i))); % number of paths for near-field

    Rh=zeros(N,N);
    for s=1:5000
        [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N, Lf, Ln, d, fc,Rmin, Rmax, kappa);
        Rh=Rh+h*h';
    end
    Rh=Rh./(5000);

    error_HF_OMP_1 = 0;  
    error_HF_SGP_1 = 0;     
    error_HF_SGP_2 = 0;         
    error_HF_OMP = 0;  
    error_HF_SGP = 0;            
    error_LS = 0;  
    error_MMSE = 0;

    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        for k = 1 : K
            [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N, Lf, Ln, d, fc,Rmin, Rmax, kappa); 
            noise = sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
            y=  h+ noise ;
            H(:,k) = h;
            Y(:,k) = y;
        end

        %% the proposed hybrid-field SGP based scheme
        
        [Hhat_HF_SGP        ]=Hybrid_SGP( Y, DFT, Polarcodebook,   Lf  *sparsity ,    Ln *sparsity, u );

        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

         %% the proposed hybrid-field SGP based scheme
        
        [Hhat_HF_SGP_1  , support , g ,sup_f ,sup_n      ]=Hybrid_SGP_1( Y, DFT, Polarcodebook,   L *sparsity_s , us ,   L);
        error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;       

        [A, G       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 10, 1, 1e-8 );
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H -  A*G   ,'fro')^2/norm( H   ,'fro')^2  ;


        %% the LS
        Hhat_LS   =  h + sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
%         Hhat_LS   =   Y  ;
        error_LS    =  error_LS +     norm(    H -  Hhat_LS  ,  'fro'   )^2   /    norm( H  ,'fro')^2  ;

        %% the hybrid-field OMP based scheme
        Hhat_homp = Hybrid_OMP( Y , DFT, Polarcodebook,    Lf*sparsity_f ,   Ln*sparsity_f  );
        error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;


        %% the hybrid-field OMP 1 1 based scheme
        Hhat_homp_1 = Hybrid_OMP_1( Y , DFT     , Polarcodebook    ,    L *sparsity_f_1  ,   L);
        error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;

       %% the MMSE
       hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2*eye(N)  )  )  )*Y; 
       error_MMSE =  error_MMSE    +     norm(    H -  hhat_MMSE  ,  'fro')^2/norm( H  ,'fro')^2  ;      


    end


NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
NMSE_LS(i) = error_LS / nbrOfSetups ;
NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;
NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;
NMSE_HF_SGP_2(i) = error_HF_SGP_2 /  nbrOfSetups ;

end

NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;
NMSE_HF_SGP_2 = 10*log10(NMSE_HF_SGP_2)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;



figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  gamma,   NMSE_LS               ,      '     k   -  o  '    , 'linewidth',  1);
plot( gamma,   NMSE_MMSE               ,      '     k   -  ^  '    , 'linewidth',  1);
plot(  gamma,   NMSE_HF_OMP               ,      '     b   -  s  '    , 'linewidth',  1);
plot(  gamma,   NMSE_HF_OMP_1               ,      '     b  --  o  '    , 'linewidth',  1);
plot(  gamma,   NMSE_HF_SGP   ,      '     r  -   p '    ,  'linewidth' , 1);
plot(  gamma,   NMSE_HF_SGP_1   ,      '     r  --   d '    ,  'linewidth' , 1);
plot(  gamma,   NMSE_HF_SGP_2   ,      '     k  --   d '    ,  'linewidth' , 2);
xlabel('$\gamma$','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
% legend(   'MMSE',  'Hybrid-field OMP (with $\gamma$) [16]' ,   'Hybrid-field OMP (without $\gamma$) [18]',     'Hybrid-field SGP (with $\gamma$)'  ,  'Hybrid-field SGP (without $\gamma$)'  ,   'Interpreter','Latex');
% ylim([ -7  6  ])

% figure('color',[1,1,1]); hold on; box on; grid on;
% % plot(  gamma,   NMSE_LS               ,      '     k   -  o  '    , 'linewidth',  1);
% plot( gamma(1:2:7),   NMSE_MMSE(1:2:7)               ,      '     k   -  ^  '    , 'linewidth',  1);
% plot(  gamma(1:2:7),   NMSE_HF_OMP(1:2:7)               ,      '     r   -  s  '    , 'linewidth',  1);
% plot(  gamma(1:2:7),   NMSE_HF_OMP_11111(1:2:7)               ,      '     r  --  >  '    , 'linewidth',  1);
% plot(  gamma(1:2:7),   NMSE_HF_SGP(1:2:7)   ,      '     g  -   p '    ,  'linewidth' , 1);
% plot(  gamma(1:2:7),   NMSE_HF_SGP_1(1:2:7)   ,      '     g  --   o '    ,  'linewidth' , 1);
% xlabel('$\gamma$','Interpreter','Latex');
% ylabel('NMSE (dB)','Interpreter','Latex');
% legend(  'MMSE',  'Hybrid-field OMP' ,   'Improved hybrid-field OMP' , 'Hybrid-field-SGP'  ,    'Improved hybrid-field SGP'  ,  'Interpreter','Latex');
% ylim([ -4  3  ])

sum(NMSE_MMSE)
sum(NMSE_HF_SGP_1)
sum(NMSE_HF_OMP_1)
