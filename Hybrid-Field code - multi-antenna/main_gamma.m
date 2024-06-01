clc;
clear
% close all
tau = 1 ;
SNR_dB =    -6  ;
SNR = SNR_dB;
SNR_linear=10.^(SNR/10.);
sigma2=1/SNR_linear;

nbrOfSetups = 1000;

%%% system parameters
N1 =    2 ; % N1 number of antenna at transmitter
N2 = 256; % N2 number of antenna at receiver

p = ones([N1,1])/(sqrt(N1));
P = diag(p);

K = 1;
L = 10; % number of all paths

sparsity_f = 8 ;
sparsity_f_1 = 4 ;

u=            0.4;       
sparsity =       1;       

% us=           0.07;   
% sparsity_s =  0.4  ;
us=           0.07;   
sparsity_s =  0.4  ;

gamma  =  0: (2/L)  :1 ;
len = length(gamma);

kappa = 10;

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space


DFT2  = (1/sqrt(N2))*exp(-1i*pi*[0:N2-1]'*[-(N2-1)/2:1:(N2/2)]*(2/N2)) ;
label_theta_R = [-(N2-1)/2:1:(N2/2)]*(2/N2);
% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5; 
[W2, label2,~, ~] = QuaCode(N2, d, lambda_c, eta, Rmin, Rmax);
S2 = size(W2, 2);

DFT1  = (1/sqrt(N1))*exp(-1i*pi*[0:N1-1]'*[-(N1-1)/2:1:(N1/2)]*(2/N1)) ;
label_theta_T = [-(N1-1)/2:1:(N1/2)]*(2/N1);
[W1, label1,~, ~] = QuaCode(N1, d, lambda_c, eta, Rmin, Rmax);
S1 = size(W1, 2);
DFT = kron( P.' * conj(DFT1) , DFT2 );
W = kron( P.' * conj(W1) , W2 );


NMSE_HF_SGP_2 = zeros(1,len );
NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );

NMSE_HF_OMP_matrix = zeros(1,len);
NMSE_HF_SGP_matrix = zeros(1,len );

NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

t0 = clock;

for  i = 1 : len

    Lf = L*gamma(i); % number of paths for far-field
    Ln =floor( L*(1-gamma(i))); % number of paths for near-field

    Rh=zeros(N2,N2);
    for s=1:5000
        [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        Rh=Rh+H*H';
    end
    Rh=Rh./(5000);

    error_LS = 0;  
    error_MMSE = 0;

    error_HF_OMP = 0;  
    error_HF_SGP = 0;     

    error_HF_OMP_1 = 0;  
    error_HF_SGP_1 = 0; 
    error_HF_SGP_2 = 0;     


    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        noise = sqrt(sigma2)*(randn(N2,N1)+1i*randn(N2,N1))/sqrt(2);
        Y=  sqrt(tau) * H *P + noise ;

        y = reshape(Y, [N1*N2,1]);

        %% the LS
%         Hhat_LS   =  H + sqrt(sigma2)*(randn(N2,N1(i))+1i*randn(N2,N1(i)))/sqrt(2);
        Hhat_LS   =   Y  ;
        error_LS    =  error_LS +     norm( H - Hhat_LS , 'fro' )^2/norm( H,'fro')^2  ;

       %% the MMSE
       hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2*eye(N2)  )  )  ) * Y; 
       error_MMSE =  error_MMSE    +     norm( H -hhat_MMSE,'fro')^2/norm( H  ,'fro')^2  ;     

        %% the hybrid-field OMP based scheme
        [Hhat_homp  ] = Hybrid_OMP( y , DFT, W,    Lf*sparsity_f ,   Ln*sparsity_f  );
        Hhat_homp = reshape(Hhat_homp, [N2,N1]);
        error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;

        %% the hybrid-field OMP 1 1 based scheme
        Hhat_homp_1 = Hybrid_OMP_1( y , DFT     , W    ,    L *sparsity_f_1  ,  L);
        Hhat_homp_1 = reshape(Hhat_homp_1, [N2,N1]);        
        error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;


        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP        ]=Hybrid_SGP( y, DFT, W,   Lf  *sparsity ,    Ln *sparsity, u );
        Hhat_HF_SGP = reshape(Hhat_HF_SGP, [N2,N1]); 
        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP_1    , support , g ,sup_f ,sup_n     ]=Hybrid_SGP_1( y, DFT, W,   L *sparsity_s , us ,   L);
        Hhat_HF_SGP_1 = reshape(Hhat_HF_SGP_1, [N2,N1]); 
        error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;

        [A, G  ]=  HF_SIGW_multiantenna(N1,  N2, y,  DFT1,  DFT2,   W1,   W2, g, support ,   sup_f ,sup_n , label1, label2, label_theta_R,  label_theta_T , fc, d, 1, 1, 0.01, 1e-8 );
        Hhat_HF_SGP_2 = A*G;
        Hhat_HF_SGP_2 = reshape(Hhat_HF_SGP_2, [N2,N1]); 
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H - Hhat_HF_SGP_2   ,'fro')^2/norm( H   ,'fro')^2  ;

    end


NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
NMSE_LS(i) = error_LS / nbrOfSetups ;
NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;
NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;
NMSE_HF_SGP_2(i) = error_HF_SGP_2 /  nbrOfSetups ;

end




NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;
NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;
NMSE_HF_SGP_2 = 10*log10(NMSE_HF_SGP_2)   ;


figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  gamma,   NMSE_LS               ,      '     k   -  o  '    , 'linewidth',  1);
plot( gamma,   NMSE_MMSE               ,      '     k   -  ^  '    , 'linewidth',  1.5);
plot(  gamma,   NMSE_HF_OMP               ,      '     b   -  s  '    , 'linewidth',   1.5);
plot(  gamma,   NMSE_HF_OMP_1               ,      '     b  --  o  '    , 'linewidth',   1.5);
plot(  gamma,   NMSE_HF_SGP   ,      '     r  -   p '    ,  'linewidth' ,  1.5);
% plot(  gamma,   NMSE_HF_SGP_1   ,      '     r  --   d '    ,  'linewidth' , 1);
plot(  gamma,   NMSE_HF_SGP_2   ,      '     r  --   d '    ,  'linewidth' ,  1.5);
xlabel('$\gamma$','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
legend(  'MMSE',  'Hybrid-field OMP (with $\gamma$) [23]' ,   'Hybrid-field OMP (without $\gamma$) [27]',     'On-grid hybrid-field SGP (with $\gamma$)'  ,  'Off-grid hybrid-field SGP (without $\gamma$)'  ,   'Interpreter','Latex');
ylim([-5 7]);



