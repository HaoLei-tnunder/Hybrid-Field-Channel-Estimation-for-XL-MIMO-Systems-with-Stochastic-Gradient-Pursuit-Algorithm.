function [hK,x,y,z,alpha] = generate_channel(N1,N2,K,P)

Xmax=P(1); Xmin=P(2); Ymax=P(3); Ymin=P(4); Zmax=P(5); Zmin=P(6);
%Xmax=0.2*P(1); Xmin=0.2*P(2); Ymax=0.2*P(3); Ymin=0.2*P(4); Zmax=0.2*P(5); Zmin=0.2*P(6);

N=N1*N2;
d=0.5;
hK=zeros(N,K);
aK=zeros(N,K);

%%%% generate the channel h_{r,k} for each user k
for k=1:K
    x = Xmax-rand*(Xmax-Xmin);
    y = Ymax-rand*(Ymax-Ymin);
    z = Zmax-rand*(Zmax-Zmin);

    alpha = (normrnd(1, .1) + 1i*normrnd(1, .1)) / sqrt(2);
    a = zeros(N,1);
    for n1=1:N1
        for n2=1:N2
            a((n1-1)*N2+n2)=alpha*exp(-1j*2*pi*sqrt((x-(n1-1-(N1-1)/2)*d)^2+(z-(n2-1-(N2-1)/2)*d)^2+y^2));
        end
    end
    hr = 1*a;
    hK(:,k) = hr./sqrt(N);
    aK(:,k) = a;
end

