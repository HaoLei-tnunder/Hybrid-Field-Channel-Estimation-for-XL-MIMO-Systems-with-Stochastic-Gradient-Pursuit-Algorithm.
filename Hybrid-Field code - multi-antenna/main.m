clc;
clear
% close all

SNR_dB = -10;%  -10  : 2  : 6  ;
len = length(SNR_dB);

nbrOfSetups = 100;

%%% system parameters
N1 = 8; % N1 number of antenna at transmitter
N2 = 256 ; % N2 number of antenna at receiver

tau = 1;
K = 1;
L = 10; % number of all paths
kappa = 10;

sparsity_f = 6 ;
sparsity_f_1 = 4 ;

% % 256   L=10   LoS = far
% %            -10       -8        -6        -4       -2          0        2         4          6         
% u=         [0.2 ,    0.4,      0.6,      0.9,      1.0,        1.8,     1.6,      1.7,        1.8  ]  ;
% sparsity = [ 1,        1,       1,        1,        1,          1,       2,         3,         4   ]  ;
% 
% %            -10       -8        -6        -4       -2          0        2         4          6    
% us=        [ 0.2 ,    0.3,      0.5,      0.8,     1.2,        1.6,     1.6,       1.9,      1.9   ]  ;
% sparsity_s =[1,        1,       1,        1,        1,         1,       1,       1,          2   ]  ;


gamma  =  6/10 ;
Lf = L*gamma; % number of paths for far-field
Ln = L -Lf; % number of paths for near-field
% M = 256; % number of pilot overhead

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space

p = ones([N1,1])/(sqrt(N1));
P = diag(p);
P1 = eye(N1)/(sqrt(N1));

DFT1  = (1/sqrt(N1))*exp(-1i*pi*[0:N1-1]'*[-(N1-1)/2:1:(N1/2)]*(2/N1)) ;
DFT2  = (1/sqrt(N2))*exp(-1i*pi*[0:N2-1]'*[-(N2-1)/2:1:(N2/2)]*(2/N2)) ;
DFT = kron( P.' * DFT1' , DFT2 );

DFT3 = DFT1';

% l = 2000;
% a = DFT(:,l) ;
% b = kron( 1/ (sqrt(N1)) * DFT3(:, floor(l/N2)+1 ) , DFT2(:,  mod(l, N2)  )  );

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5; 
[W1, ~,~, ~] = QuaCode(N1, d, lambda_c, eta, Rmin, Rmax);
[W2, ~,~, ~] = QuaCode(N2, d, lambda_c, eta, Rmin, Rmax);
S1 = size(W1, 2);
S2 = size(W2, 2);
W = kron( P.' * W1' , W2 );



Rh=zeros(N2,N2);
for s=1:5000
    [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
    Rh=Rh+H*H';
end
Rh=Rh./(5000);


NMSE_HF_OMP_1 = zeros(1,len);
NMSE_HF_SGP_1 = zeros(1,len );
NMSE_HF_OMP = zeros(1,len);
NMSE_HF_SGP = zeros(1,len );

NMSE_HF_OMP_matrix = zeros(1,len);
NMSE_HF_SGP_matrix = zeros(1,len );

NMSE_LS = zeros(1,len);
NMSE_MMSE = zeros(1,len);

NMSE_HF_SGP_2 = zeros(1,len );
NMSE_HF_SGP_3 = zeros(1,len );
NMSE_HF_SGP_4 = zeros(1,len );
NMSE_HF_SGP_5 = zeros(1,len );
NMSE_HF_SGP_6 = zeros(1,len );


t0 = clock;


for  i = 1 : len

    error_LS = 0;  
    error_MMSE = 0;

    error_HF_OMP = 0;  
    error_HF_SGP = 0;     

    error_HF_OMP_1 = 0;  
    error_HF_SGP_1 = 0; 

    error_HF_SGP_2 = 0; 
    error_HF_SGP_3 = 0; 
    error_HF_SGP_4 = 0; 
    error_HF_SGP_5 = 0; 
    error_HF_SGP_6 = 0; 


    for n = 1: nbrOfSetups

        fprintf('  iteration:[%d/%d] |  SNR:[%d/%d]  | run %.4f s\n', n,  nbrOfSetups , i ,  len , etime(clock, t0));

        SNR = SNR_dB(i);
        SNR_linear=10.^(SNR/10.);
        sigma2=1/SNR_linear;

%        [H , Hf , Hn] = generate_hybrid_field_channel_m_a(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax);
         [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1, N2 , Lf, Ln, d, fc,Rmin, Rmax, kappa);
        noise = sqrt(sigma2)*(randn(N2,N1)+1i*randn(N2,N1))/sqrt(2);
        Y=  sqrt(tau) * H * P + noise ;

        y = reshape(Y, [N1*N2,1]);

        %% the LS
%         Hhat_LS   =  H + sqrt(sigma2)*(randn(N2,N1)+1i*randn(N2,N1))/sqrt(2);
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
%         [Hhat_HF_SGP        ]=Hybrid_SGP( y, DFT, W,   Lf  *sparsity(i) ,    Ln *sparsity(i), u(i) );
%         Hhat_HF_SGP = reshape(Hhat_HF_SGP, [N2,N1]); 
%         error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;

  
        %% the proposed hybrid-field SGP based scheme
        [Hhat_HF_SGP       ]=Hybrid_SGP_1( y, DFT, W,   L  *sparsity_s(i),  us(i),   L );
        Hhat_HF_SGP = reshape(Hhat_HF_SGP, [N2,N1]); 
        error_HF_SGP =  error_HF_SGP +   norm(    H -  Hhat_HF_SGP   ,'fro')^2/norm( H   ,'fro')^2  ;
       

    end


NMSE_MMSE(i) = error_MMSE / nbrOfSetups ;
NMSE_LS(i) = error_LS / nbrOfSetups ;
NMSE_HF_OMP(i) = error_HF_OMP /   nbrOfSetups ;
NMSE_HF_SGP(i) = error_HF_SGP /  nbrOfSetups ;
NMSE_HF_OMP_1(i) = error_HF_OMP_1 /   nbrOfSetups ;


NMSE_HF_SGP_1(i) = error_HF_SGP_1 /  nbrOfSetups ;
NMSE_HF_SGP_2(i) = error_HF_SGP_2 /  nbrOfSetups ;
NMSE_HF_SGP_3(i) = error_HF_SGP_3 /  nbrOfSetups ;
NMSE_HF_SGP_4(i) = error_HF_SGP_4 /  nbrOfSetups ;
NMSE_HF_SGP_5(i) = error_HF_SGP_5 /  nbrOfSetups ;
NMSE_HF_SGP_6(i) = error_HF_SGP_6 /  nbrOfSetups ;

end



NMSE_LS = 10*log10(NMSE_LS)   ;
NMSE_MMSE = 10*log10(NMSE_MMSE)   ;
NMSE_HF_OMP = 10*log10(NMSE_HF_OMP)   ;
NMSE_HF_OMP_1 = 10*log10(NMSE_HF_OMP_1)   ;
NMSE_HF_SGP = 10*log10(NMSE_HF_SGP)   ;
NMSE_HF_SGP_1 = 10*log10(NMSE_HF_SGP_1)   ;


NMSE_HF_SGP_2 = 10*log10(NMSE_HF_SGP_2)   ;
NMSE_HF_SGP_3 = 10*log10(NMSE_HF_SGP_3)   ;
NMSE_HF_SGP_4 = 10*log10(NMSE_HF_SGP_4)   ;
NMSE_HF_SGP_5 = 10*log10(NMSE_HF_SGP_5)   ;
NMSE_HF_SGP_6 = 10*log10(NMSE_HF_SGP_6)   ;




figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  SNR_dB,   NMSE_MMSE               ,      '     k   -  o  '    , 'linewidth',  2);
% plot(  SNR_dB,   NMSE_HF_OMP               ,      '     k   --  S  '    , 'linewidth',  2);
plot(  SNR_dB,   NMSE_HF_SGP   ,      '     g  -   p '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_1   ,      '     r  -   o '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_2    ,      '     b  -   s '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_3   ,      '     c -   ^ '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_4   ,      '     r  --  > '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_5   ,      '     b  --  + '    ,  'linewidth' , 1);
plot(  SNR_dB,   NMSE_HF_SGP_6   ,      '     c  --  P '    ,  'linewidth' , 1);
xlabel('SNR (dB)','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
legend('hybrid-field-SGP' ,'hybrid-field-SGP-1'  ,   'hybrid-field-SGP-2'  , 'hybrid-field-SGP-3'  , 'hybrid-field-SGP-4'  ,'hybrid-field-SGP-5'  ,'hybrid-field-SGP-6'  ,  'Interpreter','Latex');
% legend('MMSE', 'hybrid-field-OMP' , 'hybrid-field-SGP' ,'hybrid-field-SGP-1'  ,   'hybrid-field-SGP-2'  , 'hybrid-field-SGP-3'  , 'hybrid-field-SGP-4'  ,'hybrid-field-SGP-5'  ,'hybrid-field-SGP-6'  ,  'Interpreter','Latex');

