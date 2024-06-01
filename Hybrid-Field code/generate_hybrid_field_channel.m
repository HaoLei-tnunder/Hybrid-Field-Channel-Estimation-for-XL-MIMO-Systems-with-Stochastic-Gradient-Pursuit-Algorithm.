function [H] = generate_hybrid_field_channel( N, K, L, d, fc, Rmin, Rmax, Rayleigh)

H = zeros(N,K);

for k = 1 : K

    h = zeros(N,1);


    r = rand(L,1) * (Rmax - Rmin) + Rmin; % the distance
    alpha = (randn(L,1) + 1j*randn(L,1))/sqrt(2); % the gain
    sector = pi ;
    theta = rand(L,1) * sector - sector/2; % the angle [-pi/2, pi/2]

    for l = 1 : L

        if  r(l) < Rayleigh                  % generate the near-field path components

            an = near_field_manifold( N, d, fc, r(l), theta(l));
            hn =  alpha(l) * an;
            h = h + hn ;

        else

            af = far_field_manifold(N,theta(l));
            hf =  alpha(l)*af;
            h = h + hf ;

        end


    end

    H(:,k) = h*sqrt(N/(L));

end

end









