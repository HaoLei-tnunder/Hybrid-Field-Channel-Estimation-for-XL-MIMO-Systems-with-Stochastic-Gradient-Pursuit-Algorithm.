clc;
clear
% close all

SNR_dB =    -10  ;
SNR = SNR_dB;
SNR_linear=10.^(SNR/10.);
sigma2=1/SNR_linear;

nbrOfSetups = 10;

%%% system parameters
N1 =  2 ; % N1 number of antenna at transmitter
N2 = 100:100:500; % N2 number of antenna at receiver

len = length(N2);

K = 1;
L = 10; % number of all paths

sparsity_f = 6 ;
sparsity_f_1 = 4 ;

u=            0.1;       
sparsity =       1;       

us=           0.5;   
sparsity_s =  0.2  ;


gamma  =  6/10 ;
Lf = L*gamma; % number of paths for far-field
Ln = L -Lf; % number of paths for near-field
kappa = 10;

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5; 
DFT1  = (1/sqrt(N1))*exp(-1i*pi*[0:N1-1]'*[-(N1-1)/2:1:(N1/2)]*(2/N1)) ;
[W1, ~,~, ~] = QuaCode(N1, d, lambda_c, eta, Rmin, Rmax);
S1 = size(W1, 2);




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

    p = ones([N1,1])/(sqrt(N1));
    P = diag(p);


    Rh=zeros(N2(i),N2(i));
    for s=1:5000
%         [H , Hf , Hn] = generate_hybrid_field_channel_multi_antanna(N1, N2(i) , Lf, Ln, d, fc,Rmin, Rmax);
        [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2(i) , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        Rh=Rh+H*H';
    end
    Rh=Rh./(5000);


    DFT2  = (1/sqrt(N2(i)))*exp(-1i*pi*[0:N2(i)-1]'*[-(N2(i)-1)/2:1:(N2(i)/2)]*(2/N2(i))) ;
    [W2, ~,~, ~] = QuaCode(N2(i), d, lambda_c, eta, Rmin, Rmax);
    S2 = size(W2, 2);



    DFT = kron( P.' * DFT1.' , DFT2 );
    W = kron( P.' * W1.' , W2 );
%     DFT = kron( P.' * conj(DFT1) , DFT2 );
%     W = kron( P.' * conj(W1) , W2 );


    error_LS = 0;  
    error_MMSE = 0;

    error_HF_OMP = 0;  
    error_HF_SGP = 0;     

    error_HF_OMP_1 = 0;  
    error_HF_SGP_1 = 0; 
    error_HF_OMP_matrix = 0;  
    error_HF_SGP_matrix = 0; 

    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2(i) , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        noise = sqrt(sigma2)*(randn(N2(i),N1)+1i*randn(N2(i),N1))/sqrt(2);
        Y =  H*P + noise ;

        y = reshape(Y, [N1*N2(i),1]);

        %% the LS
%         Hhat_LS   =  H + sqrt(sigma2)*(randn(N2(i),N1)+1i*randn(N2(i),N1))/sqrt(2);
        Hhat_LS   =   Y  ;
        error_LS    =  error_LS +     norm( H - Hhat_LS , 'fro' )^2/norm( H,'fro')^2  ;

       %% the MMSE
       hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2*eye(N2(i))  )  )  ) * Y; 
       error_MMSE =  error_MMSE    +     norm( H -hhat_MMSE,'fro')^2/norm( H  ,'fro')^2  ;     

        %% the hybrid-field OMP based scheme
        [Hhat_homp  ] = Hybrid_OMP( y , DFT, W,    Lf*sparsity_f ,   Ln*sparsity_f  );
        Hhat_homp = reshape(Hhat_homp, [N2(i),N1]);
        error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;

        %% the hybrid-field OMP 1 1 based scheme
        Hhat_homp_1 = Hybrid_OMP_1( y , DFT     , W    ,    L *sparsity_f_1  ,  L);
        Hhat_homp_1 = reshape(Hhat_homp_1, [N2(i),N1]);        
        error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;


        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP        ]=Hybrid_SGP( y, DFT, W,   Lf  *sparsity ,    Ln *sparsity, u );
        Hhat_HF_SGP = reshape(Hhat_HF_SGP, [N2(i),N1]); 
        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP_1       ]=Hybrid_SGP_1( y, DFT, W,   L *sparsity_s , us ,   L);
        Hhat_HF_SGP_1 = reshape(Hhat_HF_SGP_1, [N2(i),N1]); 
        error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;
       
    

    end


NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
NMSE_LS(i) = error_LS / nbrOfSetups ;
NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;
NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;


end




NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;
NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;



figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  N1,   NMSE_LS               ,      '     k   -  +  '    , 'linewidth',  1);
plot(  N2,   NMSE_MMSE               ,      '       k   --  ^  '    , 'linewidth',  2);
plot(  N2,   NMSE_HF_OMP               ,      '     k   -  s  '    , 'linewidth',  1);
plot(  N2,   NMSE_HF_OMP_1               ,      '     B   -  >  '    , 'linewidth',  1);
plot(  N2,   NMSE_HF_SGP   ,      '     R  -   p '    ,  'linewidth' , 1);
plot(  N2,   NMSE_HF_SGP_1   ,      '     g  -   o '    ,  'linewidth' , 1);
xlabel('Number of antennas at the UE','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
legend(  'Hybrid-field OMP (with $\gamma$) [16]',  'Hybrid-field OMP (without $\gamma$) [18]','Hybrid-field SGP (with $\gamma$)',  'Hybrid-field SGP (without $\gamma$)', 'Interpreter','Latex');


