function [H,H_LOS,H_NLOS,r_n1,r_n2,theta_n1,theta_n2] = generate_near_field_MIMO_channel(N1, N2, theta, phi, r, lambda,Ln, Rmin, Rmax,beta)
%N1 number of antenna at transmitter
%N2 number of antenna at receiver
%beta = 3;
% generate LOS path
c = 3.0e8;
% H_LOS = zeros(N2,N1);
H_NLOS = zeros(N2,N1);
% R = zeros(N2,N1);
H_LOS = generate_LOS_near_field_MIMO_channel(N1, N2, theta,phi,r,lambda)';

% H_LOS = H_LOS./norm(H_LOS);
    
% generate NLOS path 
d = lambda/2;
fc = lambda/c;
r_n1 = rand(Ln,1) * (Rmax - Rmin) + Rmin;
r_n2 = rand(Ln,1) * (Rmax - Rmin) + Rmin;% the distance 
sector = 2*pi/3; 
theta_n1 = rand(Ln,1) * sector - sector/2; % the angle [-pi/3, pi/3]
alpha_n1 = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain 

% alpha_n1 =1./(r_n1.*r_n2) .* alpha_n1;

theta_n2 = rand(Ln,1) * sector - sector/2; % the angle [-pi/3, pi/3]

% alpha_n2 = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain 
for l = 1:Ln
   an1 = near_field_manifold( N1, d, fc, r_n1(l), theta_n1(l));
   an2 = near_field_manifold( N2, d, fc, r_n2(l), theta_n2(l));
   h = an2 * an1';
%    H_NLOS = H_NLOS + alpha_n1(l)*alpha_n2(l) * h;
    H_NLOS = H_NLOS + alpha_n1(l) * h;
end 
% H_NLOS = H_NLOS/sqrt(Ln);
H = H_LOS+H_NLOS;
% H_NLOS = H_NLOS./norm(H_NLOS);
% H = sqrt(beta/(1+beta))*H_LOS+sqrt(1/(1+beta))*H_NLOS;



% % generate the near-field path components
% r_n = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance 
% sector = 2*pi/3; 
% theta_n = rand(Ln,1) * sector - sector/2; % the angle [-pi/3, pi/3]
% alpha_n = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain 
% for l = 1:Ln
%    an = near_field_manifold( N, d, fc, r_n(l), theta_n(l));
%    h = h + alpha_n(l) * an;
% end 
% 
% % generate the far-field path components
% alpha_f = (randn(Lf,1) + 1j*randn(Lf,1))/sqrt(2); % the gain 
% sector = 2*pi/3; 
% theta_f = rand(Lf,1) * sector - sector/2; % the angle [-pi/3, pi/3]
% for l = 1:Lf
%     af = far_field_manifold(N,theta_f(l));
%     h = h + alpha_f(l)*af;
% end
% 
% h = h*sqrt(N/(Ln+Lf));

