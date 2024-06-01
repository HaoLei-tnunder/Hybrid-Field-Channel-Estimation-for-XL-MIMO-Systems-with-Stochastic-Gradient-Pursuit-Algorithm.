function[hhat]=Hybrid_OMP_1(y,Uf,Un,  t , L )

% Af =P* Uf;
% An =P* Un;
Af = Uf;
An = Un;

D = [  Af, An  ] ;

[M,N]=size(Af);
r=y ;
R=zeros(M,  t  )         ; 

posf=[];
for i=1:t
    R(:,i) = r ;
	Pf = Af'*r;
    [~,pos]=max(abs(Pf));
	posf=[posf,pos(end)];
	hhat_f=zeros(N,1);
	hhat_f(posf)=pinv(Af(:,posf))*y;
    r=y-Af*hhat_f;  
end

posn=[];
[~,S]=size(An);
rm = r;
hhat_m = zeros(N+S,1);
hhat_m(1:N) = hhat_f;
for gamma = ((L-1)/L) : (-1/L) :0
    posn  =  posf( 1:    floor(gamma*t)    );
%     posn  =  posf( 1:    floor(gamma*t)+1    );    
    r = R( :  , max(1, floor(gamma*t + 1  ))   );   % ceil
    for i  =  1  :   (  floor(  t *(1-gamma))   )
	    Pn = An'*r;
        [~,pos]=max(abs(Pn));
% 	    posn=[posn,pos(end)];
	    posn=[posn,pos(end)+N];        
	    hhat_1=zeros(  N+S  ,1);
	    hhat_1(posn)=pinv(D(:,posn))*y;
        r=y - D*hhat_1;  
    end
    if norm(r)<norm(rm)
        rm = r;
        hhat_m = hhat_1;
        sss = posn;
    end
end


hhat = Uf*hhat_m(1:N)  +Un*hhat_m(N+1:N+S);


D = [  Af  ,  An  ] ;

support =  find(hhat_m);
G =hhat_m(support);
h1 = D(:,support ) * G;
sup_f = find(hhat_f);

sup_n = find(hhat_m(N+1:N+S)) - 256;

a = 1;

    
        