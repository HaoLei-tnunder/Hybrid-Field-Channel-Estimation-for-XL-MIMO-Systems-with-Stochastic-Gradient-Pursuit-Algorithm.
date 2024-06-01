function [H] = generate_hybrid_field_channel_4(N, Lf, Ln, d, fc, kappa)

H = zeros(N*N,1);
L = Lf+ Ln;
c = 3e8;

x = 500-rand(L,1)*(500-10);
y = 500-rand(L,1)*(500-10);
z = 500-rand(L,1)*(500-10);

for i = 1 : L
    alpha = (randn(1,1) + 1j*randn(1,1))/sqrt(2); % the gain 
    a = zeros(N*N,1);
    for n1=1:N
        for n2=1:N
            a((n1-1)*N+n2)=alpha*exp(-1j*2*pi*fc/c*sqrt((x(i)-(n1-1-(N-1)/2)*d)^2+(z(i)-(n2-1-(N-1)/2)*d)^2+y(i)^2));
        end
    end
    H = H + sqrt(1/(1+kappa))* a ./sqrt(N*N/L) ;

end

x = 500-rand(1,1)*(500-500);
y = 500-rand(1,1)*(500-500);
z = 500-rand(1,1)*(500-500);

alpha = (randn(1,1) + 1j*randn(1,1))/sqrt(2); % the gain 
a = zeros(N*N,1);
for n1=1:N
    for n2=1:N
        a((n1-1)*N+n2)=alpha*exp(-1j*2*pi*fc/c*sqrt((x(i)-(n1-1-(N-1)/2)*d)^2+(z(i)-(n2-1-(N-1)/2)*d)^2+y(i)^2));
    end
end
H = H + sqrt(kappa/(1+kappa))* a ./sqrt(N*N) ;