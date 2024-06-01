function [   Hhat  , support , g ,sup_f ,sup_n]=Hybrid_SGP_1( Y, DFT, Polarcodebook, t ,  u ,   L)


% Af = P* DFT;
% An = P*Polarcodebook;
Af = DFT;
An = Polarcodebook;

[~,N]=size(Af);
[M,S]=size(An);

%%% initialization
r = Y;
support_F = [];
support_N = [];
hhat_f  =  zeros(N,1);
R=zeros( M  , t )         ;
D = [  Af  ,  An  ] ;

for i = 1 : t
    R(:,i) = r  ;
    Product_F =   Af'  *  r   ;
    [~,posf]=   max(abs(Product_F))  ;
    support_F =  [support_F ,posf(end)]   ;
    x = hhat_f(support_F)  ;
    for j = 1 : M
        B = Af(:,support_F) ;
        c = B(j,:);
        dj = Y(j) ;
        ej = dj -  c * x  ;
        x = x + u * ej * B(j,:)' ;
    end
    hhat_f(support_F) =  x ;
    r = Y - Af*hhat_f;
end

rm = r;
hhat_m = zeros(N+S,1);
hhat_n = zeros(N+S,1);
hhat_m(1:N) = hhat_f;


for gamma = ((L-1)/L) : (-1/L) :0
    support_N  =  support_F( 1:    floor(gamma*t)    );
    r = R( :  ,  floor(gamma*t + 1  ) );   % ceil
    for i  =  1  :   (  floor(  t*(1-gamma))   )
        Product_N =An'  *  r   ;
        [~,posn]=   max(abs(Product_N))  ;
%         support_N =  [support_N ,posn(end)]   ;
        support_N =  [support_N ,posn(end) + N]   ;            
        xx = hhat_n(support_N)  ;
        for j = 1 : M
            B = D(:,support_N) ;
            c = B(j,:);
            dj = Y(j) ;
            ej = dj -  c * xx  ;
            xx = xx + u * ej * B(j,:)' ;
        end
        hhat_n(support_N) =  xx ;
        r = Y - D*hhat_n    ;
    end
    if norm(r)<norm(rm)
        rm = r;
        hhat_m = hhat_n;
    end

end

Hhat = DFT*hhat_m(1:N)  +Polarcodebook*hhat_m(N+1:N+S);

support = find(hhat_m);
g = hhat_m(support);


hhat_far =  hhat_m(1:N);
hhat_near = hhat_m(N+1:N+S);

sup_f = find(hhat_far);
g_f = hhat_far(sup_f );

sup_n = find(hhat_near);
g_n = hhat_near(sup_n );

