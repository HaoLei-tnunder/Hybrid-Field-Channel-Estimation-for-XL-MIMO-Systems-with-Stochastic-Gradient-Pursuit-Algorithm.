function[ hhat, support] = SGP( y , A , u   ,  L )
%y:measurement 
%A:sensing matrix
%t:number of iteration

[M,Q]=size(A);        

r=y;  %Initialization of residual

support=[];
hhat=zeros(Q,1);


for l = 1 :  L

	Product = A'*r;                                 %   N*Q                                                                      ;                N*N
    [~,pos]=max(abs(Product));                %   Q                                                   ;                    N
	support=[support,pos(end)];            %       1                                        ;                1

    x = hhat(support)  ;

% 	hhat=zeros(N,1);

    for j = 1 : M

        B = A(:,support) ;       

        c = B(j,:);  

        dj = y(j) ;              

        ej = dj -  c * x  ;                   %            L              ;                                L

        x = x + u * ej * B(j,:)' ;                %            L              ;                                L

    end

	hhat(support) =  x ;

    r = y - A * hhat;  

end

