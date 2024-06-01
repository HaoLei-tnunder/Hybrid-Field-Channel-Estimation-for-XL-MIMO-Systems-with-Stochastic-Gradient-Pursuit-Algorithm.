function [h,hf,hn] = generate_hybrid_field_channel_3(N, Lf, Ln, d, fc,Rmin, Rmax, kappa)

% generate the far-field path components
hf = zeros(N,1);

alpha_f = (randn(Lf,1) + 1j*randn(Lf,1))/sqrt(2); % the gain
sector = 3*pi/3;
theta_f = rand(Lf,1) * sector - sector/2; % the angle [-pi/2, pi/2]
r_f = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance
for l = 1:Lf
    rou = rand(1,1) * (4 - 2) + 2; 
    beta_f = 0.001*((r_f(l))^(-rou));    
    af = far_field_manifold(N,theta_f(l));
%     hf = hf +sqrt(beta_f) * sqrt(1/(1+kappa))*alpha_f(l)*af;
    hf = hf +  sqrt(1/(1+kappa))*alpha_f(l)*af;    
end

% generate the near-field path components
hn = zeros(N,1);
r_n = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance
sector = 3*pi/3;
theta_n = rand(Ln,1) * sector - sector/2; % the angle [-pi/2, pi/2]
alpha_n = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain
for l = 1:Ln
    rou = rand(1,1) * (4 - 2) + 2; 
    beta_n = 0.001*((r_n(l))^(-rou));
    an = near_field_manifold( N, d, fc, r_n(l), theta_n(l));
%     hn = hn + sqrt(beta_n) * sqrt(1/(1+kappa))*alpha_n(l) * an;
    hn = hn +  sqrt(1/(1+kappa))*alpha_n(l) * an;    
end


h = hn+hf;
h = h*sqrt(N/(Ln+Lf));

% %% LoS = far
% sector = 3*pi/3;
% theta = rand(1,1) * sector - sector/2; % the angle [-pi/2, pi/2]
% a = far_field_manifold(N,theta);
% h = h + sqrt(N)*sqrt(kappa/(1+kappa))*a;



%% LoS = near
r = rand(1,1) * (Rmax - Rmin) + Rmin; % the distance
sector = 3*pi/3;
rou = rand(1,1) * (4 - 2) + 2;
beta = 0.001*((r)^(-rou));
theta = rand(1,1) * sector - sector/2; % the angle [-pi/2, pi/2]
a = near_field_manifold( N, d, fc, r, theta);
% h = h + sqrt(beta) *  sqrt(N)*sqrt(kappa/(1+kappa))*a;
h = h +   sqrt(N)*sqrt(kappa/(1+kappa))*a;
% h=  beta * h;

