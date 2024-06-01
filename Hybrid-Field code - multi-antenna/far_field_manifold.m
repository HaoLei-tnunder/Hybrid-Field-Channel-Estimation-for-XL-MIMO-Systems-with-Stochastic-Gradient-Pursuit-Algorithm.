function  at = far_field_manifold(Nt,theta)
    at = exp(-1i*pi*[0:Nt-1]'*sin(theta));
    at = at / norm(at);
end

