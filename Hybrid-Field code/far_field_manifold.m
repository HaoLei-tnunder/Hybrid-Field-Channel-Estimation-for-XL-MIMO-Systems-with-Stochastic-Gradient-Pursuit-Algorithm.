function  at = far_field_manifold(Nt,theta)
    at = exp(-1i*pi*[0:Nt-1]'*sin(theta)) ;
%     at1 = exp(-1i*pi*[-(Nt - 1)/2:(Nt - 1)/2]'*sin(theta));    
%     at = exp(-1i*pi*[0:Nt-1]'*sin(theta))./ sqrt(Nt) ;
    at = at / norm(at);
end

