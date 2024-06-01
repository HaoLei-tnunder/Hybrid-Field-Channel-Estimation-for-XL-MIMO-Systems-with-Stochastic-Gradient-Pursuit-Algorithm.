function [dict, label, dict_cell, label_cell] = QuaCode(N, d, lambda, eta, rho_min, rho_max)
    c = 3e8;
    theta = -1 + 2/N : 2/N : 1;
    dict = cell(N, 1);
    label = cell(N, 1);
    Zmax = (N * d)^2 / 2 / lambda / eta^2;
    kmax = floor(Zmax/rho_min);
    for idx = 1:N
        Z = (N * d)^2 * ( 1 - theta(idx)^2) / 2 / lambda / eta^2;
        kmax = floor(Z/rho_min);
        kmin = floor(Z/rho_max) + 1;
        
        r = zeros(1, kmax - kmin + 2);
        r(:,1) = (N * d)^2 * 2 / lambda;
        r(:,2:end) = Z./(kmin:kmax);
        
        dict{idx} = zeros(N, kmax + 1);
        label{idx} = zeros(2, kmax + 1);
%         disp('-------')
%         disp(kmax + 1)
%         disp(kmax - kmin + 2)
        for t = 1 : kmax - kmin + 2
            dict{idx}(:, t) = near_field_manifold( N, d, c/lambda, r(t), asin(theta(idx)) );
            label{idx}(:, t) = [theta(idx), r(t)]';
        end
    end
    dict_cell = dict;
    label_cell = label;
    dict = merge(dict, N, N);
    label = merge(label, N, 2);
end

function B = merge(A, N, Q)
    S = zeros(1, N);
    for idx = 1:N
        S(idx) = size(A{idx}, 2);
    end
    B = zeros(Q, sum(S));
    for idx = 1:N
        B(:, sum(S(1:idx)) - S(idx) + 1: sum(S(1:idx))) = A{idx};
    end
end