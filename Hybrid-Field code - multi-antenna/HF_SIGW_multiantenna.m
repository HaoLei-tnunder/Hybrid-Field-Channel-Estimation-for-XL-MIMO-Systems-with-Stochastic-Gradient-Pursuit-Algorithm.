function [A, G  ] = HF_SIGW_multiantenna(N1,N2,Y, DFT1,  DFT2,   W1,   W2, g0, support ,   sup_f ,sup_n , label1, label2, label_theta_R,  label_theta_T, f, d,  N_iter, lr_alpha, lr_th, th)
    p = ones([N1,1])/(sqrt(N1));
    P = diag(p);
% [N, ~] = size(Y);

S2 = size(W2, 2);

L = length(support);
Lf = length(sup_f);
Ln = length(sup_n);
theta_fr = zeros(1, Lf);
theta_ft = zeros(1, Lf);
sup_fr = zeros(Lf,1);
sup_ft = zeros(Lf,1);

for i = 1 : Lf
    a = floor( sup_f(i) / N2  ) +1;
    if a == N1+1
        a = N1;
    end     
    b = mod(   sup_f(i),   N2   );
    if b == 0
        b = N2;
    end
    theta_fr(i) = label_theta_R(b);
    theta_ft(i) = label_theta_T(a);    
    sup_fr(i) = b;
    sup_ft(i) = a;    
end

theta_nr = zeros(1, Ln);
theta_nt = zeros(1, Ln);
r_nr = zeros(1, Ln);
r_nt = zeros(1, Ln);
sup_nr = zeros(Ln,1);
sup_nt = zeros(Ln,1);

for i = 1 : Ln
    a = floor( sup_n(i) / S2  ) +1;
    if a == N1+1
        a = N1;
    end    
    b = mod(   sup_n(i),   S2   );
    if b == 0
        b = S2;
    end
    theta_nr(i) = label2(1,b);
    theta_nt(i) = label1(1,a);    
    r_nr(i) = label2(2,b);
    r_nt(i) = label1(2,a);    
    sup_nr(i) = b;
    sup_nt(i) = a;        
end

DFT11 = DFT1(:,   sup_ft);
DFT22 = DFT2(:,   sup_fr);

W11 = W1(:,  sup_nt);
W22 = W2(:,  sup_nr);
A = zeros(N1*N2,  L);

for i = 1 : Lf
    A(:, i) = kron( sqrt(1/N1) * conj(DFT11(:, i)) , DFT22(:, i) );
end
for i = 1 + Lf :  L
    A(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , W22(:, i-Lf) );
end

DFT = kron( P.' * conj(DFT1) , DFT2 );
W = kron( P.' * conj(W1) , W2 );
Af = DFT;
An = W;
D = [  Af  ,  An  ] ;
A0 = D(:, support);

% Ini

G =g0;
G_prev = G;

alpha_r = (1 - theta_nr.^2) ./ 2 ./ r_nr;
alpha_t = (1 - theta_nt.^2) ./ 2 ./ r_nt;
c = 3e8;


% lambda = 10;

Res = Y - A0 * g0 ;
Res_list = [];
Rnorm = norm(Res, 'fro')^2;
Res_list = [Res_list, Rnorm];

c1 = 0.1;
cl_alpha = 0.01;
% calculate gradient


for iter = 1:N_iter
    %% update theta -------  theta_fr
    dtheta01 = df_dtheta_fr(Y,  A,   Lf,theta_fr, DFT11,N1, N2   );
    direction_theta = dtheta01;
    lr_th_iter = lr_th;
    A_iter = A;
    while lr_th_iter > 1e-6
        theta_iter = theta_fr - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
        rowf = (0:N2-1)';
        % update A
        A_iter_f = exp(- 1j*  pi * rowf * theta_iter )/sqrt(N2);
        for i = 1 : Lf
            A_iter(:, i) = kron( sqrt(1/N1) * conj(DFT11(:, i)) , A_iter_f (:, i) );
        end
        G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y; 
        Res2 = Y - A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_th_iter * norm(direction_theta, 'fro')^2
            break;
        else
            lr_th_iter = lr_th_iter * 0.8;
        end
     end

    % update theta(A) and gains(G)
    theta_fr = theta_fr - lr_th_iter * direction_theta;
    theta_fr(theta_fr > 1) = theta_fr(theta_fr > 1) - 2;
    theta_fr(theta_fr < -1) = theta_fr(theta_fr < -1) + 2;
    % update A
    A_iter_f = exp(- 1j*  pi * rowf * theta_fr )/sqrt(N2);
    for i = 1 : Lf
        A(:, i) = kron( sqrt(1/N1) * conj(DFT11(:, i)) , A_iter_f (:, i) );
    end

    %% update theta -------  theta_ft
    dtheta02 = df_dtheta_ft(Y,  A,   Lf,theta_ft, DFT22,N1 );
    direction_theta = dtheta02;
    lr_th_iter = lr_th;
    A_iter = A;
    while lr_th_iter > 1e-6
        theta_iter = theta_ft - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
        rowf = (0:N1-1)';
        % update A
        A_iter_f = exp(- 1j*  pi * rowf * theta_iter )/sqrt(N1);
        for i = 1 : Lf
            A_iter(:, i) = kron( sqrt(1/N1) * conj(A_iter_f (:, i)) , DFT22(:, i) );
        end
        G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y; 
        Res2 = Y - A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_th_iter * norm(direction_theta, 'fro')^2
            break;
        else
            lr_th_iter = lr_th_iter * 0.8;
        end
     end

    % update theta(A) and gains(G)
    theta_ft = theta_ft - lr_th_iter * direction_theta;
    theta_ft(theta_ft > 1) = theta_ft(theta_ft > 1) - 2;
    theta_ft(theta_ft < -1) = theta_ft(theta_ft < -1) + 2;
    % update A
    A_iter_f = exp(- 1j*  pi * rowf * theta_ft )/sqrt(N1);
    for i = 1 : Lf
        A_iter(:, i) = kron( sqrt(1/N1) * conj(A_iter_f (:, i)) , DFT22(:, i) );
    end


    %% update theta -------  theta_nr
    dtheta03 = df_dtheta_nr(Y,  A,   Lf,theta_nr, W11,N1, N2   , alpha_r, d, f );
    direction_theta = dtheta03;
    lr_th_iter = lr_th;
    A_iter = A;
    row = (-(N2 - 1)/2:(N2 - 1)/2)' ;
    while lr_th_iter > 1e-6
        theta_iter = theta_nr - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
                      
        % update A     
        rn_iter =   - ( row * d ) .* theta_iter + ( row * d ).^2 .* alpha_r;         
        A_iter_n = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N2);
        for i = 1 + Lf :  L
            A_iter(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_iter_n(:, i-Lf) );
        end
        G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y; 
        Res2 = Y - A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_th_iter * norm(direction_theta, 'fro')^2
            break;
        else
            lr_th_iter = lr_th_iter * 0.8;
        end
     end

    % update theta(A) and gains(G)
    theta_nr = theta_nr - lr_th_iter * direction_theta;
    theta_nr(theta_nr > 1) = theta_nr(theta_nr > 1) - 2;
    theta_nr(theta_nr < -1) = theta_nr(theta_nr < -1) + 2;

    % update A
    rn_iter =   - ( row * d ) .* theta_nr + ( row * d ).^2 .* alpha_r;         
    A_iter_n = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N2);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_iter_n(:, i-Lf) );
    end


    %% update theta -------  theta_nt
    dtheta04 = df_dtheta_nt(Y,  A,   Lf,  theta_nt, W22,N1,  alpha_t, d, f );
    direction_theta = dtheta04;
    lr_th_iter = lr_th;
    A_iter = A;
    row = (-(N1 - 1)/2:(N1 - 1)/2)' ;
    while lr_th_iter > 1e-6
        theta_iter =  theta_nt - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
                      
        % update A     
        rn_iter =   - ( row * d ) .* theta_iter + ( row * d ).^2 .* alpha_t;         
        A_iter_n = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N1);
        for i = 1 + Lf :  L
            A_iter(:, i) = kron( sqrt(1/N1) * conj(A_iter_n(:, i-Lf)) , W22(:, i-Lf) );
        end
        G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y; 
        Res2 = Y - A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_th_iter * norm(direction_theta, 'fro')^2
            break;
        else
            lr_th_iter = lr_th_iter * 0.8;
        end
     end

    % update theta(A) and gains(G)
    theta_nt = theta_nt - lr_th_iter * direction_theta;
    theta_nt(theta_nt > 1) = theta_nt(theta_nt > 1) - 2;
    theta_nt(theta_nt < -1) = theta_nt(theta_nt < -1) + 2;

    % update A
    rn_iter =   - ( row * d ) .* theta_nt + ( row * d ).^2 .* alpha_t;         
    A_iter_n = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N1);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(A_iter_n(:, i-Lf)) , W22(:, i-Lf) );
    end

    %% update distance  alpha  ----------------   r_nr

    if L > Lf
        dalpha0 = df_dalpha_r(Y , A, Lf , theta_nr, W11, N1,N2, alpha_r , d, f );
        direction_alpha = dalpha0;
        lr_alpha_iter = lr_alpha;
        A_iter = A;  
        row = (-(N2 - 1)/2:(N2 - 1)/2)' ;        
        while lr_alpha_iter > 1e-6  
            alpha_iter = alpha_r - lr_alpha_iter * direction_alpha;
            alpha_iter( alpha_iter < 0 ) = 1e-10;
    
            % update A
            rn_iter =  - ( row * d ) .* theta_nr + ( row * d ).^2 .* alpha_iter;
            A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N2);            
            for i = 1 + Lf :  L
                A_iter(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_n_iter(:, i-Lf) );
            end
    
            G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y;
            
            Res2 = Y - A_iter*G_iter;
            Rnorm2 = norm(Res2, 'fro')^2;
            if Rnorm2 < Rnorm - cl_alpha * lr_alpha_iter * norm(direction_alpha, 'fro')^2
                break;
            else
                lr_alpha_iter = lr_alpha_iter * 0.8;
            end
        end
         
        % update distanace(r) and gains(G)
        alpha_r  = alpha_r  - lr_alpha_iter * direction_alpha;
        alpha_r ( alpha_r  <= 0 ) = 1e-10;
    
        % update A
            rn_iter =  - ( row * d ) .* theta_nr + ( row * d ).^2 .* alpha_iter;
            A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N2);            
            for i = 1 + Lf :  L
                A(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_n_iter(:, i-Lf) );
            end
           
    end

    %% update distance  alpha  ----------------   r_nt

    if L > Lf
        dalpha1 = df_dalpha_t(Y , A, Lf , theta_nt, W22, N1,alpha_t , d, f );
        direction_alpha = dalpha1;
        lr_alpha_iter = lr_alpha;
        A_iter = A; 
        row = (-(N1 - 1)/2:(N1 - 1)/2)' ;                
        while lr_alpha_iter > 1e-6  
            alpha_iter = alpha_t - lr_alpha_iter * direction_alpha;
            alpha_iter( alpha_iter < 0 ) = 1e-10;
    
            % update A
            rn_iter =  - ( row * d ) .* theta_nt + ( row * d ).^2 .* alpha_iter;
            A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N1);            
            for i = 1 + Lf :  L
                A_iter(:, i) = kron( sqrt(1/N1) * conj(A_n_iter(:, i-Lf)) , W22(:, i-Lf) );
            end
            G_iter = pinv( A_iter'*A_iter ) * A_iter' * Y;    
            Res2 = Y - A_iter*G_iter;
            Rnorm2 = norm(Res2, 'fro')^2;
            if Rnorm2 < Rnorm - cl_alpha * lr_alpha_iter * norm(direction_alpha, 'fro')^2
                break;
            else
                lr_alpha_iter = lr_alpha_iter * 0.8;
            end
        end
         
        % update distanace(r) and gains(G)
        alpha_r  = alpha_r  - lr_alpha_iter * direction_alpha;
        alpha_r ( alpha_r  <= 0 ) = 1e-10;
    
        % update A
            rn_iter =  - ( row * d ) .* theta_nt + ( row * d ).^2 .* alpha_r;
            A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N1);            
            for i = 1 + Lf :  L
                A_iter(:, i) = kron( sqrt(1/N1) * conj(A_n_iter(:, i-Lf)) , W22(:, i-Lf) );
            end
       
        R =  A;
        G = pinv( R'*R  ) * R' * Y;
    
        % obtain res and res norm
        Res = Y - A*G;
        Rnorm = norm(Res, 'fro')^2;
        Res_list = [Res_list, Rnorm];
        
    end

    % update epsilon
    gamma = norm(G - G_prev, 'fro');
  
    % early stopping
    if gamma < th
       break; 
    end
    % update previous G
    G_prev = G;  
end

end


%%

function dff = df_dtheta_fr(Y,  A,   Lf,theta_fr, DFT11,N1, N2   )
    dA = dA_dtheta_fr(A, Lf , theta_fr, DFT11, N1,N2);
    R =   A;
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(1:Lf);
end

function dA = dA_dtheta_fr(A, Lf , theta_fr, DFT11, N1,N2)
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (0:N2-1)' ;    
    A_iter_f = (-1j*pi*row).* exp(- 1j*  pi * row * theta_fr )/sqrt(N2);
    for i = 1 : Lf
        dA(:, i) = kron( sqrt(1/N1) * conj(DFT11(:, i)) , A_iter_f (:, i) );
    end    
end

%%

function dff = df_dtheta_ft(Y,  A,   Lf,theta_ft, DFT22,N1 )
    dA = dA_dtheta_ft(A, Lf , theta_ft, DFT22, N1);
    R =   A;
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(1:Lf);    
end

function dA = dA_dtheta_ft(A, Lf , theta_ft, DFT22, N1)
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (0:N1-1)' ;    
    A_iter_f = (-1j*pi*row).* exp(- 1j*  pi * row * theta_ft )/sqrt(N1);
    for i = 1 : Lf
        dA(:, i) = kron( sqrt(1/N1) * conj(A_iter_f(:, i)) , DFT22(:, i) );
    end    
end

%%

function dff = df_dtheta_nr(Y,  A,   Lf,theta_nr, W11,N1, N2   , alpha_r, d, f )
    dA = dA_dtheta_nr(A, Lf , theta_nr, W11, N1,N2, alpha_r , d, f );
    R =   A;
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(Lf+1:L);
end

function dA = dA_dtheta_nr(A, Lf , theta_nr, W11, N1,N2, alpha_r, d, f )
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (-(N2- 1)/2:(N2 - 1)/2)' ;
     c = 3e8;
    rn_iter =   - ( row * d ) .* theta_nr + ( row * d ).^2 .* alpha_r;         
    A_iter_n =  (1j*2*pi*f/c).*( row * d )  .* exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N2);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_iter_n(:, i-Lf) );
    end

end

%%

function dff = df_dtheta_nt(Y,  A,   Lf,  theta_nt, W22,N1,  alpha_t, d, f )
    dA = dA_dtheta_nt(A, Lf , theta_nt, W22, N1,alpha_t , d, f );
    R =   A;
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l );
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(Lf+1:L);    
end

function dA = dA_dtheta_nt(A, Lf , theta_nt, W22, N1,alpha_t, d, f  )
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (-(N1- 1)/2:(N1 - 1)/2)' ;
     c = 3e8;
    rn_iter =   - ( row * d ) .* theta_nt + ( row * d ).^2 .* alpha_t;         
    A_iter_n =  (1j*2*pi*f/c).*( row * d )  .* exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N1);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(A_iter_n(:, i-Lf)) , W22(:, i-Lf) );
    end

end

%%

function dff = df_dalpha_r(Y , A, Lf , theta_nr, W11, N1,N2, alpha_r , d, f )
    dA = dA_dalpha_r(A, Lf , theta_nr, W11, N1,N2, alpha_r, d , f );
    R = A;
    E = R'*R;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l );
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(Lf+1:L);
end


function dA = dA_dalpha_r(A, Lf , theta_nr, W11, N1,N2, alpha_r, d, f  )
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (-(N2- 1)/2:(N2 - 1)/2)' ;
     c = 3e8;
    rn_iter =   - ( row * d ) .* theta_nr + ( row * d ).^2 .* alpha_r;         
    A_iter_n =  (1j*2*pi*f/c).*( - (row * d).^2  ) .* exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N2);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(W11(:, i-Lf)) , A_iter_n(:, i-Lf) );
    end
end

%%

function dff = df_dalpha_t(Y , A, Lf , theta_nt, W22, N1,alpha_t , d, f )
    dA = dA_dalpha_t(A, Lf , theta_nt, W22, N1, alpha_t , d, f );
    R = A;
    E = R'*R;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l );
        dR_dl = dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
    dff = df(Lf+1:L);
end


function dA = dA_dalpha_t(A, Lf , theta_nt, W22, N1, alpha_t, d , f )
    [N, L] = size(A);
    dA = zeros(N, L);
    row = (-(N1- 1)/2:(N1 - 1)/2)' ;
     c = 3e8;
    rn_iter =   - ( row * d ) .* theta_nt + ( row * d ).^2 .* alpha_t;         
    A_iter_n =  (1j*2*pi*f/c).*( - (row * d).^2  ) .* exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N1);
    for i = 1 + Lf :  L
        A(:, i) = kron( sqrt(1/N1) * conj(A_iter_n(:, i-Lf)) , W22(:, i-Lf) );
    end
end






