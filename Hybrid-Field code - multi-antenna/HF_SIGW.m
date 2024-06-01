function [A, G  ] = HF_SIGW(Y, DFT, Polarcodebook, g0, support, r_n0, theta_n0,  theta_f0, f, d, L_min, pruning, N_iter, lr_alpha, lr_th, th)

% [N, ~] = size(Y);

Af = DFT;
An = Polarcodebook;
D = [  Af  ,  An  ] ;
A0 = D(:, support);

[N, L] = size(A0);

Lf = length(theta_f0);
% Ln = L - Lf ;

% Ini
A = A0;

G =g0;
G_prev = G;

theta_f = theta_f0;
theta_n = theta_n0;
theta = [theta_f ,theta_n ];
r = r_n0;
alpha = (1 - theta_n.^2) ./ 2 ./ r;

c = 3e8;
row = (-(N - 1)/2:(N - 1)/2)' ;
rowf = (0:N-1)' ;

% lambda = 10;

Res = Y - A0 * g0 ;
Res_list = [];
Rnorm = norm(Res, 'fro')^2;
Res_list = [Res_list, Rnorm];

c1 = 0.1;
cl_alpha = 0.01;
% calculate gradient


for iter = 1:N_iter
    %% update  theta
    dtheta0 = df_dtheta(Y,  A,  r ,  Lf,f,d);
    direction_theta = dtheta0;
    lr_th_iter = lr_th;
    
    while lr_th_iter > 1e-6  
        theta_iter = theta - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
        
                
        % update A

        rn_iter =   - ( row * d ) .* theta_iter( Lf+1:L ) + ( row * d ).^2 .* alpha;  
%         A_iter_f =  exp( 1j*  pi * row * theta_iter( 1:Lf ) )/sqrt(N);
        A_iter_f = exp(- 1j*  pi * rowf * theta_iter( 1:Lf ) )/sqrt(N);
        A_iter_n = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N);        
        A_iter = [A_iter_f, A_iter_n];
%         for i = 1 : L
%             A_iter(:,i) = A_iter(:,i) / norm(A_iter(:,i));
%         end

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
    theta = theta - lr_th_iter * direction_theta;
    theta(theta > 1) = theta(theta > 1) - 2;
    theta(theta < -1) = theta(theta < -1) + 2;

    
    % update A
    rn =  - ( row * d ) .* theta( Lf+1:L ) + ( row * d ).^2 .* alpha;
%     A_f =  exp( 1j*  pi * row * theta_iter( 1:Lf ) )/sqrt(N);
    A_f = exp(- 1j*  pi * rowf * theta( 1:Lf ) )/sqrt(N);
    A_n = exp( - 1j*2*pi*f*( rn )/c ) / sqrt(N);
    A = [A_f, A_n];
%     for i = 1 : L
%             A(:,i) = A(:,i) / norm(A(:,i));
%     end
    G = pinv( A'*A ) * A' * Y;
    
    %% update distance  alpha

    if L > Lf
        dalpha0 = df_dalpha(Y, A,f,d, Lf, L);
        direction_alpha = dalpha0;
        lr_alpha_iter = lr_alpha;
        while lr_alpha_iter > 1e-6  
            alpha_iter = alpha - lr_alpha_iter * direction_alpha;
            alpha_iter( alpha_iter < 0 ) = 1e-10;
    
            % update A
            rn_iter =  - ( row * d ) .* theta( Lf+1:L ) + ( row * d ).^2 .* alpha_iter;
%             A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c);
            A_n_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N);            
            A_iter = [A_f,A_n_iter ];
%             for i = 1 : L
%                 A_iter(:,i) = A_iter(:,i) / norm(A_iter(:,i));
%             end    
    
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
        alpha = alpha - lr_alpha_iter * direction_alpha;
        alpha( alpha <= 0 ) = 1e-10;
    
        % update A
        rn =  - ( row * d ) .* theta( Lf+1:L ) + ( row * d ).^2 .* alpha;
%         A_n = exp( - 1j*2*pi*f*( rn )/c) ;
        A_n = exp( - 1j*2*pi*f*( rn )/c) / sqrt(N);        
        A = [A_f, A_n];    
%         for i = 1 : L
%             A(:,i) = A(:,i) / norm(A(:,i));
%         end    

        R =  A;
        G = pinv( R'*R  ) * R' * Y;
    
        % obtain res and res norm
        Res = Y - A*G;
        Rnorm = norm(Res, 'fro')^2;
        Res_list = [Res_list, Rnorm];
        
    end

    % update epsilon
    gamma = norm(G - G_prev, 'fro');
  
%     % pruning
% %     if L > L_min && iter > 10
%     if L > L_min
%         Gnorm = sum(abs(G).^2, 2);
%         index = find(Gnorm > pruning);
% %         index = find( Gnorm > pruning * max(Gnorm));
%         r = r(index);
%         theta = theta(index);
%         alpha = alpha(index);
%         A = A(:, index);
%         G = G(index, :);
%         dtheta0 = dtheta0(index);
%         dalpha0 = dalpha0(index);
%         direction_theta = direction_theta(index);
%         direction_alpha = direction_alpha(index);
%         L = numel(index);
%     end


    % early stopping
    if gamma < th
       break; 
    end
    % update previous G
    G_prev = G;  
end

end

function df = df_dtheta(Y,  A,  r ,  Lf,f,d)
    dA = dA_dtheta(A, Lf , r,f,d);
    R =   A;
%     size(R'*R)
%     size(D)
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
end


function dff = df_dalpha(Y,  A, f,d, Lf, L)
% f = - tr{ Y'R(R'R + D/lambda)^(-1)R'Y}
% R = W'A
    dA = dA_dalpha(A,  f,d, Lf, L);
    R = A;
    E = R'*R;
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


function dA = dA_dtheta(A, Lf ,   r,f,d)

    [N, L] = size(A);
    dA = zeros(N, L);
%     row = (-(N - 1)/2:(N - 1)/2)' ;
    row = (0:N-1)' ;    
  
    dA(:, 1:Lf) = (-1j*pi*row).*A(:, 1:Lf);

    c = 3e8;
%     row = (-(N - 1)/2:(N - 1)/2)' ;
%     rn = sqrt( r.^2 + ( row * d ).^2 - 2 * row * d * ( r.*theta ));    

    dA(:,  Lf+1 : L) =(1j*2*pi*f/c).*( row * d )  .* A(:  ,  Lf+1 : L) ;         

end


function dA = dA_dalpha(A, f,d, Lf, L)
    N = size(A, 1);
    dA = zeros(N, L);
    c = 3e8;
    row = (-(N - 1)/2:(N - 1)/2)' ;
%     dA = (-1j*2*pi*f/c).*( r./ rn - row * d * theta ./rn - 1) .* A;
    dA(:, 1:Lf) = 0;
    dA(:,  Lf+1 : L) = (1j*2*pi*f/c).*( - (row * d).^2  ) .* A(:,  Lf+1 : L) ;              
end


