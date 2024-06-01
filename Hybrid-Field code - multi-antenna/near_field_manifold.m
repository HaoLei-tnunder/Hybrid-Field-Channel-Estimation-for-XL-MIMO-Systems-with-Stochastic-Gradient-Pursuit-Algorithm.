function  at = near_field_manifold( Nt, d, f, r0, theta0 )
    c = 3e8;
    lambda=2*d;
    nn = [-(Nt-1)/2:1:(Nt-1)/2]';
     r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sin(theta0));
%     at = exp(-1j*2*pi*f*(r - r0)/c) * r0./r;
%     at = at / norm(at);
     a = (-1j*2*pi*f*(r - r0)/c);
    b = -1j* (  nn.*pi*sin(theta0)- (nn).^2*pi/4/r0/lambda*((cos(theta0))^2)*d^2  ) ;
    at = exp(a) ;
%     at = exp(a)* r0./ r ;    
    at = at / norm(at);


end
