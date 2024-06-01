function[  SE  ]  =  functionComputeMonteCarloSE( Hhat , H , Pu , sigma2_u)


[Nb,Nu]=size(H);

%V =(     (   Hhat * Pu*(Pu') * (Hhat') + sigma2_u *   eye(Nb) )^(-1)   ) * Hhat * Pu;                 %MMSE
V = Hhat ;   %MR

%n =  sqrt(sigma2_u)*(randn(Nb,1)+1i*randn(Nb,1))/sqrt(2);


D = (V') * H*Pu ;

%Term1 = (V') * n * (n') * V;
Term1 =  sigma2_u *  (V')  * V;

Term2 =    (D')  * (Term1^(-1))    *  D  ;

Term3 =   eye(Nu) +   Term2    ;

Term4  =  log2 (  det(      Term3     ) ) ;

SE    =    real(    Term4      );



