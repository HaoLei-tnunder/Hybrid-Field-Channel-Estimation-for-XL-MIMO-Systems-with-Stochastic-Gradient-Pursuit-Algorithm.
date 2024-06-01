function [   Hhat , support , g ,sup_f ,sup_n]=Hybrid_SGP( Y, DFT, Polarcodebook, Lf   ,  Ln , u )

% Af = P* DFT;
% An = P*Polarcodebook;
Af = DFT;
An = Polarcodebook;

[~,N]=size(Af);   
[M,S]=size(An);   

%%% initialization
R = Y;
support_F = [];
support_N = [];
hhat_f  =  zeros(N,1);
hhat_n   =  zeros(S,1);

a = 5;
b = 0.03;
cc = 3;



for i = 1 : Lf
	Product_F =Af'  *  R   ;                           %  N*N                               
    [~,posf]=   max(abs(Product_F))  ;                 %  
	support_F =  [support_F ,posf(end)]   ;            %
    x = hhat_f(support_F)  ;                           %
    for j = 1 : M        %  N(2L+1)
        B = Af(:,support_F) ;                          %
        c = B(j,:);                                    %
        dj = Y(j) ;                                    %
        ej = dj -  c * x  ;                            %  L

%         u = b*(  3/(   1  + exp(  ( -a * abs( atan( ej  ) )^cc      )        )         )    -1   );
        x = x + u * ej * B(j,:)' ;                     %  L
    end
	hhat_f(support_F) =  x ;                           %
    R = Y - Af*hhat_f;                                 %  N*N
end



for i = 1 :   Ln
	Product_N =An'  *  R   ;                           %  SN                                            
    [~,posn]=   max(abs(Product_N))  ;                 % 
	support_N =  [support_N ,posn(end)]   ;            %
    xx = hhat_n(support_N)  ;                          %
    for j = 1 : M                       %  N(2L+1)
        B = An(:,support_N) ;                          %    
        c = B(j,:);                                    %
        dj = Y(j) ;                                    %
        ej = dj -  c * xx  ;                           %
%         u = b*(  3/(   1  + exp(  ( -a * abs( atan( ej  ) )^cc      )        )         )    -1   );        
        xx = xx + u * ej * B(j,:)' ;                   %
    end
	hhat_n(support_N) =  xx ;                          %
    if length(support_F)>0                             %
        R = Y - Af*hhat_f   -      An*hhat_n    ;      %  
    else
        R = Y -   An*hhat_n    ;                       %  SN
    end
end


if length(support_N)>0
    if length(support_F)>0
        Hhat = DFT*hhat_f   +  Polarcodebook*hhat_n;   %  N*N+SN
    else
        Hhat =   Polarcodebook*hhat_n;
    end
else
    Hhat = DFT*hhat_f   ;
end


D = [  Af  ,  An  ] ;
hhat =  transpose([     transpose(hhat_f), transpose(hhat_n) ]);
support =  find(hhat);
g =hhat(support);

sup_f = find(hhat_f);

sup_n = find(hhat_n);

% H = D(:,support) * hhat(support);
% H1 = D * hhat;
% a = 1;


