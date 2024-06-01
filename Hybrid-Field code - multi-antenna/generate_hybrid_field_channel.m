function [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax)

% generate the far-field path components
hf = zeros(N,1);
alpha_f = (randn(Lf,1) + 1j*randn(Lf,1))/sqrt(2); % the gain 
sector = 3*pi/3; 
theta_f = rand(Lf,1) * sector - sector/2; % the angle [-pi/2, pi/2]
for l = 1:Lf
    af = far_field_manifold(N,theta_f(l));
    hf = hf + alpha_f(l)*af;
end

% generate the near-field path components
hn = zeros(N,1);
r_n = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance 
sector = 3*pi/3; 
theta_n = rand(Ln,1) * sector - sector/2; % the angle [-pi/2, pi/2]
alpha_n = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain 
for l = 1:Ln
   an = near_field_manifold( N, d, fc, r_n(l), theta_n(l));
   hn = hn + alpha_n(l) * an;
end 


h = hn+hf;
h = h*sqrt(N/(Ln+Lf));

