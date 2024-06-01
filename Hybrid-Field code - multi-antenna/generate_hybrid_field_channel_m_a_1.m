function [H , Hf , Hn] = generate_hybrid_field_channel_m_a_1(N1,N2, Lf, Ln, d, fc,Rmin, Rmax,kappa )


% generate the far-field path components
Hf = zeros(N2,N1);
alpha_f = (randn(Lf,1) + 1j*randn(Lf,1))/sqrt(2); % the gain
sector = 3*pi/3;
theta_f1 = rand(Lf,1) * sector - sector/2; % the angle [-pi/2, pi/2]
theta_f2 = rand(Lf,1) * sector - sector/2; % the angle [-pi/2, pi/2]

for l = 1:Lf
    af1 = far_field_manifold(N1,theta_f1(l));
    af2 = far_field_manifold(N2,theta_f2(l));
    Hf = Hf + sqrt(1/(1+kappa))*alpha_f(l)*af2*af1';
end

% generate the near-field path components
Hn = zeros(N2,N1);
r_n1 = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance
r_n2 = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance
sector = 3*pi/3;
theta_n1 = rand(Ln,1) * sector - sector/2; % the angle [-pi/2, pi/2]
theta_n2 = rand(Ln,1) * sector - sector/2; % the angle [-pi/2, pi/2]
alpha_n = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain
for l = 1:Ln
    an1 = near_field_manifold( N1, d, fc, r_n1(l), theta_n1(l));
    an2 = near_field_manifold( N2, d, fc, r_n2(l), theta_n2(l));
    Hn = Hn + sqrt(1/(1+kappa))*alpha_n(l) * an2*an1';
end

H = Hn+Hf;
H = H*sqrt(N1*N2/(Ln+Lf));


sector = 3*pi/3;
theta_1 = rand(1,1) * sector - sector/2; % the angle [-pi/2, pi/2]
theta_2 = rand(1,1) * sector - sector/2; % the angle [-pi/2, pi/2]

% a1 = ones([N1,1]);
a1 = far_field_manifold(N1,theta_1);
a2 = far_field_manifold(N2,theta_2);
H = H + sqrt(N1*N2)*sqrt(kappa/(1+kappa))*a2*a1';
% H = H + sqrt(N1*N2/(Ln+Lf))*sqrt(kappa/(1+kappa))*a2*a1';


% H = H*sqrt(N/(Ln+Lf));

