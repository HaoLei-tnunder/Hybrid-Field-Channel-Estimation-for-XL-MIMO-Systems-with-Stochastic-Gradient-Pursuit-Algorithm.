function[hhat,   x]=OMP_1(y,A,t)
%y:measurement
%A:sensing matrix
%t:number of iteration

[~,N]=size(A);     %   N*Q                                               ; N*N
r=y;  %Initialization of residual
support=[];

x = 0;

while norm(r) > 60
    for i=1:t
    	Product = A'*r;                     %   N*Q                                                                      ;                N*N
        [~,pos]=max(abs(Product));                     %   Q                                                   ;                    N
    	support=[support,pos(end)];                         %       1                                        ;                1
    	hhat=zeros(N,1);

    	hhat(support)=  pinv(   A(:,support)   )  *   y;                   %   N*L*L + N*L                                ;              N*L*L + N*L

        r=y-A*hhat;                                              %                  N*Q                                               ;              N*N
    end
    x=x+1;
end

