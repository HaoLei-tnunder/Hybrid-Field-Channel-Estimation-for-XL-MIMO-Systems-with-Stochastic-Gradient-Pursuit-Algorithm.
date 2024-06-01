clc;
clear
% close all

SNR_dB =  -6  : 2  : 10  ;
len = length(SNR_dB);

nbrOfSetups = 1000;

%%% system parameters
N = 256; % number of beams (transmit antennas)
% M = 256;
% N_RF=4;
% Q = 64;  % number of pilot blocks
K = 1;
L   = 10; % number of all paths
kappa = 10;

sparsity_f = 8 ;
sparsity_f_1 = 4 ;
sparsity_f_2 = 6 ;




u=             0.4      ;
sparsity =    1      ;


%                   
us=             0.03       ;
sparsity_s =  0.6     ;


% test
% %                  -6        -4       -2          0        2            4          6         8           10
% u=         [    0.4,      0.4,      0.8,        0.8,     1.2,      1.2,        1.3,      1.3        1.3     ]  ;
% sparsity = [   1,        1,        1,          1,       3,          4,           5,        5            5      ]  ;
% test
% %                   -6           -4        -2          0           2             4          6            8         10     
% us=        [     0.03 ,    0.03,     0.06,      0.08,     0.16,        0.35,     0.9,       1.1         1.3      ]  ;
% sparsity_s =[  0.6,       0.6,       0.6,       0.4,       0.4,          0.4,       0.4,        1,           1.5         ]  ;

     
gamma  = 5/10  ;
Lf =floor( L*gamma); % number of paths for far-field
Ln = L -Lf; % number of paths for near-field
% M = 256; % number of pilot overhead

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space


% DFT1  = (1/sqrt(N))*exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N)) ;
s = 1;
D = s*N; %字典规模
% row = (-(N - 1)/2:(N - 1)/2)' ;
row =  (0:N-1)'  ;
% col = -1 + 2/D : 2/D : 1 ;
col = [-(N-1)/2:1:(N/2)]*(2/N) ;
DFT  =  exp( -1j*  pi * row * col ) / sqrt(N);

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
% Rayleigh=50;
eta = 2.5;
% [Polarcodebook, ~,~, ~] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
[Polarcodebook, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
S = size(Polarcodebook, 2);
% label_theta = -1 + 2/N : 2/N : 1;
label_theta = [-(N-1)/2:1:(N/2)]*(2/N);



% Rh=zeros(N,N);
% % Rhp=zeros(S,S);
% % Rha=zeros(N,N);
% for s=1:5000
% %     [h      ]=generate_hybrid_field_channel( N, K, L, d, fc, Rmin, Rmax, Rayleigh);    
% %     [h   ,hf,hn   ]=generate_hybrid_field_channel_1(N, Lf, Ln, d, fc,Rmin, Rmax);
%     [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N, Lf, Ln, d, fc,Rmin, Rmax, kappa);
% 
% %     h = channel_norm(h);
%     Rh=Rh+h*h';
% %     NORM_h = NORM_h + norm( h  ,'fro')^2 / 5000;
%     %     hp = x1*h;
%     %     Rhp = Rhp + hp*hp';
%     %     ha = inv(DFT)*h;
%     %     Rha = Rha + ha*ha';
% end
% Rh=Rh./(5000);



NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );

NMSE_OMP = zeros(1,len);

NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

NMSE_HF_SGP_2 = zeros(1,len );
NMSE_HF_SGP_3 = zeros(1,len );
NMSE_HF_SGP_4 = zeros(1,len );
NMSE_HF_SGP_5 = zeros(1,len );
NMSE_HF_SGP_6 = zeros(1,len );

H = zeros(N,K) ;
% Y = zeros(M,K) ;
Y = zeros(N,K) ;

t0 = clock;

T = zeros(len, 1);

for  i = 1 : len

    error_HF_OMP = 0;
    error_HF_SGP = 0;


    error_LS = 0;
    error_MMSE = 0;
    error_HF_OMP_1 = 0;
    error_HF_SGP_1 = 0;

    error_OMP = 0;

        error_HF_SGP_2 = 0;
        error_HF_SGP_3 = 0;
        error_HF_SGP_4 = 0;
        error_HF_SGP_5 = 0;
        error_HF_SGP_6 = 0;


    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        % H = generate_hybrid_field_channel( N, K, L, d, fc, Rmin, Rmax, Rayleigh);



        %         P = (2*randi(2,M,N) - 3)/sqrt(N);
        %         P = ((rand(M,N)>0.5)*2-1)/sqrt(M);


%         [h , hf , hn] = generate_hybrid_field_channel_1(N, Lf, Ln, d, fc,Rmin, Rmax);
        [h   ,hf,hn   ]=generate_hybrid_field_channel_2(N, Lf, Ln, d, fc,Rmin, Rmax, kappa);        
%         h = channel_norm(h);
        SNR = SNR_dB(i);
        SNR_linear=10.^(SNR/10.);
        sigma2= 1/SNR_linear;       
%         sigma2=  1/SNR_linear; 
        noise = sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
        %             noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        y=  h+ noise ;
        H = h;
        Y = y;



%         %% the LS
%         Hhat_LS   =  h + sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
%         %         Hhat_LS   =   Y  ;
%         error_LS    =  error_LS +     norm(    H -  Hhat_LS  ,  'fro'   )^2   /    norm( H  ,'fro')^2  ;
% 
%         %% the MMSE
%         hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2*eye(N)  )  )  ) * Hhat_LS;
%         error_MMSE =  error_MMSE    +     norm(    H -  hhat_MMSE  ,  'fro')^2/norm( H  ,'fro')^2  ;


%         %% the far-field OMP based scheme
% 
%         [Hhat_omp  ] =  OMP( Y , DFT,    L*sparsity_f_2  );
%         error_OMP =  error_OMP +     norm(    H -  Hhat_omp  ,  'fro')^2/norm( H  ,'fro')^2  ;
% 
% 
%         %% the hybrid-field OMP based scheme
% 
% 
%         [Hhat_homp  ] = Hybrid_OMP( Y , DFT, Polarcodebook,    Lf*sparsity_f ,   Ln*sparsity_f  );
%         error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;
% 
%         %% the hybrid-field OMP 1 1 based scheme
%         Hhat_homp_1 = Hybrid_OMP_1( Y , DFT     , Polarcodebook    ,    L *sparsity_f_1  ,  L);
%         error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;



        %% the proposed hybrid-field SGP based scheme
        [ Hhat_HF_SGP   , support , g ,sup_f ,sup_n   ]=Hybrid_SGP( Y, DFT, Polarcodebook,   Lf  *sparsity(i)  ,    Ln *sparsity(i) , u(i)  );        
%         [Hhat_HF_SGP        ]=Hybrid_SGP( Y, DFT, Polarcodebook,   Lf  *sparsity(i) ,    Ln *sparsity(i), u(i) );
        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

        [A1, G1       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 10, 0.1, 1e-8 );
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H -  A1*G1   ,'fro')^2/norm( H   ,'fro')^2  ;        

        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP_1   , support , g ,sup_f ,sup_n      ]=Hybrid_SGP_1( Y, DFT, Polarcodebook,   L *sparsity_s(i)  , us(i)  ,   L);       
        error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;

        [A, G       ]=  HF_SIGW( Y  ,  DFT, Polarcodebook, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 10, 1, 1e-8 );
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H -  A*G   ,'fro')^2/norm( H   ,'fro')^2  ;


  
    end


    NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
    NMSE_LS(i) = error_LS / nbrOfSetups ;
    NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
    NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
    NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;

    NMSE_OMP(i) = error_OMP /   nbrOfSetups ;

    NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;
    NMSE_HF_SGP_2(i) = error_HF_SGP_2 /  nbrOfSetups ;
    NMSE_HF_SGP_3(i) = error_HF_SGP_3 /  nbrOfSetups ;
    NMSE_HF_SGP_4(i) = error_HF_SGP_4 /  nbrOfSetups ;
    NMSE_HF_SGP_5(i) = error_HF_SGP_5 /  nbrOfSetups ;
    NMSE_HF_SGP_6(i) = error_HF_SGP_6 /  nbrOfSetups ;

end



NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;
NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;

NMSE_OMP = 10*log10(NMSE_OMP)   ;

NMSE_HF_SGP_2 = 10*log10(NMSE_HF_SGP_2)   ;
NMSE_HF_SGP_3 = 10*log10(NMSE_HF_SGP_3)   ;
NMSE_HF_SGP_4 = 10*log10(NMSE_HF_SGP_4)   ;
NMSE_HF_SGP_5 = 10*log10(NMSE_HF_SGP_5)   ;
NMSE_HF_SGP_6 = 10*log10(NMSE_HF_SGP_6)   ;


CustomCNN1 = [ 0.73485489,    -2.23060515,    -5.8820587,  -10.2209164 ,    -14.90587973,       -18.82588284,       -21.37318898,        -22.86507087,        -23.78381633    ];              %     SNR = 0


figure('color',[1,1,1]); hold on; box on; grid on;
plot(  SNR_dB,   NMSE_LS               ,      '     k   -  +  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_MMSE               ,      '       b   -  ^  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_OMP               ,      '     m   -  >  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_HF_OMP               ,      '     r   -  s  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_HF_SGP   ,      '     g  -   p '    ,  'linewidth' , 1.5);
xlabel('SNR (dB)','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
legend('LS', 'MMSE' ,'Far-field OMP [36]'   , 'Hybrid-field OMP (with $\gamma$) [23]'  ,'On-grid hybrid-field SGP (with $\gamma$)' , 'Interpreter','Latex');
% legend( 'LS', 'MMSE' , 'Hybrid-field OMP'  ,'Hybrid-field SGP'  ,     'Interpreter','Latex');


figure('color',[1,1,1]); hold on; box on; grid on;
plot(  SNR_dB,   NMSE_OMP               ,      '     m   -  ^  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_HF_OMP               ,      '     k   -  s  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_HF_OMP_1               ,      '     b   -  >  '    , 'linewidth',  1.5);
plot(  SNR_dB,   NMSE_HF_SGP   ,      '     r  -   p '    ,  'linewidth' , 1.5);
% plot(  SNR_dB,   NMSE_HF_SGP_1   ,      '     g  -   o '    ,  'linewidth' , 1.5);
plot(  SNR_dB,   NMSE_HF_SGP_2   ,      '     g  -   o '    ,  'linewidth' , 1.5);
plot(  SNR_dB,   CustomCNN1   ,      '     k  -   o '    ,  'linewidth' , 3);
xlabel('SNR (dB)','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
% legend('Far-field OMP [34]'   , 'Hybrid-field OMP (with $\gamma$) [22]'  ,'Hybrid-field OMP (without $\gamma$) [26]'  ,  'On-grid hybrid-field SGP (with $\gamma$)' ,  'Off-grid hybrid-field SGP (without $\gamma$)'  , 'Interpreter','Latex');
legend('Far-field OMP [36]'   , 'Hybrid-field OMP (with $\gamma$) [23]'  ,'Hybrid-field OMP (without $\gamma$) [27]'  ,  'On-grid hybrid-field SGP (with $\gamma$)' ,  'Off-grid hybrid-field SGP (without $\gamma$)' ,'CustomCNN', 'Interpreter','Latex');


