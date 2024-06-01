function[hhat  ]=Hybrid_OMP(y,Uf,Un,Lf,Ln )

% Af = P*Uf;
% An = P*Un;
Af =  Uf;
An = Un;

[M,N]=size(Af);

r=y;  %Initialization of residual
posf=[];
for i=1:Lf
	Pf = Af'*r;
    [~,pos]=max(abs(Pf));
	posf=[posf,pos(end)];
	hhat_f=zeros(N,1);
	hhat_f(posf)=pinv(Af(:,posf))*y;
    r=y-Af*hhat_f;  
end

posn=[];
[M,N]=size(An);
for i=1:Ln
	Pn = An'*r;
    [~,pos]=max(abs(Pn));
	posn=[posn,pos(end)];
	hhat_n=zeros(N,1);
	hhat_n(posn)=pinv(An(:,posn))*y;
    if length(posf)>0
        r=y-An*hhat_n-Af*hhat_f;  
    else
        r=y-An*hhat_n;
    end
end


if length(posn)>0 
    if length(posf)>0
        hhat = Uf*hhat_f+Un*hhat_n;
    else
        hhat = Un*hhat_n;
    end
else
    hhat = Uf*hhat_f;
end
    
        