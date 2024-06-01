clc;
clear

tau = 1 ;
SNR_dB =  -6  : 2  : 6  ;
len=length(SNR_dB);


nbrOfSetups = 1000;

%%% system parameters
N1 =    1 ; % N1 number of antenna at transmitter
N2 = 256; % N2 number of antenna at receiver

pt = ones([N1,1])/(sqrt(N1));
Pt = diag(pt);

pu =  ones([N1,1])/(sqrt(N1));
Pu = diag(pu);

% %Communication bandwidth
% B = 20e6;
%
% %Noise figure (in dB)
% noiseFigure = 7;
%
% %Compute noise power
% noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
% % sigma2 = db2pow(noiseVariancedBm-30);



K = 1;
L = 10; % number of all paths
kappa = 10;

sparsity_f = 8 ;
sparsity_f_1 = 4 ;


% %   test           
% u=         [     1,      1,        1,          1,       1.2,      1.2,        1.3,      1.3        1.3     ]  ;
% sparsity = [   1,     1,        1,          1,       3,          4,           5,        5            5      ]  ;
% 
%  
% %    test
% us=        [     0.06 ,    0.07,     0.1,      0.1,     0.16,        0.35,     0.9,       1.1         1.3      ]  ;
% sparsity_s =[ 0.4,       0.4,       0.4,       0.4,       0.4,          0.4,       0.4,        1,           1.5         ]  ;


u=          1 ;
sparsity =   1  ;
           
us=            0.1  ;
sparsity_s =0.4 ;



gamma  =  4/10 ;



fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space


DFT2  = (1/sqrt(N2))*exp(-1i*pi*[0:N2-1]'*[-(N2-1)/2:1:(N2/2)]*(2/N2)) ;
DFT1  = (1/sqrt(N1))*exp(-1i*pi*[0:N1-1]'*[-(N1-1)/2:1:(N1/2)]*(2/N1)) ;

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5;
[W2, label, dict_cell, label_cell] = QuaCode(N2, d, lambda_c, eta, Rmin, Rmax);
S2 = size(W2, 2);
[W1, ~,~, ~] = QuaCode(N1, d, lambda_c, eta, Rmin, Rmax);
S1 = size(W1, 2);

DFT = kron( Pt.' * DFT1.' , DFT2 );
W = kron( Pt.' * W1.' , W2 );
label_theta = [-(N2-1)/2:1:(N2/2)]*(2/N2);

t0 = clock;

Lf = L*gamma; % number of paths for far-field
Ln =floor( L*(1-gamma)); % number of paths for near-field

Rh=zeros(N2,N2);
for s=1:5000
    [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
    Rh=Rh+H*H';
end
Rh=Rh./(5000);

SE_Monte_MR_Combining_perfectCSI = zeros(len,1);
SE_Monte_MR_Combining_LS = zeros(len,1);
SE_Monte_MR_Combining_MMSE = zeros(len,1);
SE_Monte_MR_Combining_HF_OMP = zeros(len,1);
SE_Monte_MR_Combining_HF_OMP_1 =zeros(len,1);
SE_Monte_MR_Combining_HF_SGP = zeros(len,1);
SE_Monte_MR_Combining_HF_SGP_1 = zeros(len,1);
SE_Monte_MR_Combining_HF_SGP_2 = zeros(len,1);


NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );

NMSE_OMP = zeros(1,len);

NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

NMSE_HF_SGP_2 = zeros(1,len );
% NMSE_HF_SGP_3 = zeros(1,len );
% NMSE_HF_SGP_4 = zeros(1,len );
% NMSE_HF_SGP_5 = zeros(1,len );
% NMSE_HF_SGP_6 = zeros(1,len );



for i = 1:len

error_LS = 0;
error_MMSE = 0;
error_HF_OMP = 0;
error_HF_SGP = 0;
error_HF_OMP_1 = 0;
error_HF_SGP_1 = 0;
error_HF_SGP_2 = 0;

    SNR = SNR_dB(i);
    SNR_linear=10.^(SNR/10.);

    sigma2_t= 1 /SNR_linear;

    sigma2_u = sigma2_t;


    for n = 1: nbrOfSetups


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         channel estimation         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        noise = sqrt(sigma2_t)*(randn(N2,N1)+1i*randn(N2,N1))/sqrt(2);
        Y=  sqrt(tau) * H *Pt + noise ;

        y = reshape(Y, [N1*N2,1]);

        %---------------- the LS
        Hhat_LS   =   Y  ;
            error_LS    =  error_LS +     norm( H - Hhat_LS , 'fro' )^2/norm( H,'fro')^2  ;

        %---------------- the MMSE
        hhat_MMSE=   Rh* (  inv  (  Rh+(sigma2_t*eye(N2)  )  )  ) * Y;
            error_MMSE =  error_MMSE    +     norm( H -hhat_MMSE,'fro')^2/norm( H  ,'fro')^2  ;

        %------------------ the hybrid-field OMP based scheme
        [Hhat_homp  ] = Hybrid_OMP( y , DFT, W,    Lf*sparsity_f ,   Ln*sparsity_f  );
        Hhat_homp = reshape(Hhat_homp, [N2,N1]);
            error_HF_OMP =  error_HF_OMP +     norm(    H -  Hhat_homp  ,  'fro')^2/norm( H  ,'fro')^2  ;

        %--------------------- the hybrid-field OMP 1 1 based scheme
        Hhat_homp_1 = Hybrid_OMP_1( y , DFT     , W    ,    L *sparsity_f_1  ,  L);
        Hhat_homp_1 = reshape(Hhat_homp_1, [N2,N1]);
            error_HF_OMP_1 =  error_HF_OMP_1 +     norm(    H -  Hhat_homp_1  ,  'fro')^2/norm( H  ,'fro')^2  ;


        %------------------------ the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP        ]=Hybrid_SGP( y, DFT, W,   Lf  *sparsity(i) ,    Ln *sparsity(i), u(i) );
        Hhat_HF_SGP = reshape(Hhat_HF_SGP, [N2,N1]);
            error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

        %------------------------- the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP_1   , support , g ,sup_f ,sup_n    ]=Hybrid_SGP_1( y, DFT, W,   L *sparsity_s(i) , us(i) ,   L);
        Hhat_HF_SGP_1 = reshape(Hhat_HF_SGP_1, [N2,N1]);
            error_HF_SGP_1 =  error_HF_SGP_1 +   norm(    H -  Hhat_HF_SGP_1   ,'fro')^2/norm( H   ,'fro')^2  ;

        [A, G       ]=  HF_SIGW( Y  ,  DFT, W, g, support  , label(2, sup_n), label(1, sup_n), label_theta(1,sup_f ) , fc, d, L, 3, 1, 0.1, 0.01, 1e-8 );
        Hhat_HF_SGP_2 = A*G;
        error_HF_SGP_2 =  error_HF_SGP_2 +   norm(    H - Hhat_HF_SGP_2   ,'fro')^2/norm( H   ,'fro')^2  ;


        %%         SE


        [ SE_Monte_MR_Combining_perfectCSI(i) ]  =    SE_Monte_MR_Combining_perfectCSI(i)   +    functionComputeMonteCarloSE( H , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_LS(i) ]  =                SE_Monte_MR_Combining_LS(i)                +    functionComputeMonteCarloSE( Hhat_LS , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_MMSE(i) ]  =         SE_Monte_MR_Combining_MMSE(i)         + functionComputeMonteCarloSE( hhat_MMSE , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_HF_OMP(i) ]  =     SE_Monte_MR_Combining_HF_OMP(i)      +  functionComputeMonteCarloSE( Hhat_homp , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_HF_OMP_1(i) ]  = SE_Monte_MR_Combining_HF_OMP_1(i)   +  functionComputeMonteCarloSE( Hhat_homp_1 , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_HF_SGP(i) ]  =      SE_Monte_MR_Combining_HF_SGP(i)        + functionComputeMonteCarloSE( Hhat_HF_SGP , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_HF_SGP_1(i) ]  =  SE_Monte_MR_Combining_HF_SGP_1(i)     +  functionComputeMonteCarloSE( Hhat_HF_SGP_1 , H , Pu , sigma2_u)/  nbrOfSetups    ;
        [ SE_Monte_MR_Combining_HF_SGP_2(i) ]  =  SE_Monte_MR_Combining_HF_SGP_2(i)     +  functionComputeMonteCarloSE( Hhat_HF_SGP_2 , H , Pu , sigma2_u)/  nbrOfSetups    ;

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
plot( SNR_dB,   SE_Monte_MR_Combining_perfectCSI               ,      '     k   -  ^  '    , 'linewidth',  1);
plot(  SNR_dB,   SE_Monte_MR_Combining_HF_SGP_2               ,      '      g  --   d  '    , 'linewidth',  1);
% plot(  SNR_dB,   SE_Monte_MR_Combining_HF_SGP_1   ,      '     g  --   d '    ,  'linewidth' , 1);
plot(  SNR_dB,   SE_Monte_MR_Combining_HF_SGP   ,      '     g -   p '    ,  'linewidth' , 1);
plot(  SNR_dB,   SE_Monte_MR_Combining_HF_OMP_1   ,      '     r  --   > '    ,  'linewidth' , 1);
plot(  SNR_dB,   SE_Monte_MR_Combining_HF_OMP   ,      '     r  -   O '    ,  'linewidth' , 1);
plot(  SNR_dB,   SE_Monte_MR_Combining_LS               ,      '     b   -  s  '    , 'linewidth',  1);
% plot(  SNR_dB,   SE_Monte_MR_Combining_MMSE               ,      '     b  --  o  '    , 'linewidth',  1);
xlabel('SNR (dB)','Interpreter','Latex');
ylabel('Achievable Rate (bit/s/Hz)','Interpreter','Latex');
legend(  'Perfect CSI',  'Off-grid hybrid-field SGP (without $\gamma$)' , 'On-grid hybrid-field SGP (with $\gamma$)' ,   'Hybrid-field OMP (without $\gamma$) [27]',  'Hybrid-field OMP (with $\gamma$) [23]', 'LS','Interpreter','Latex');



