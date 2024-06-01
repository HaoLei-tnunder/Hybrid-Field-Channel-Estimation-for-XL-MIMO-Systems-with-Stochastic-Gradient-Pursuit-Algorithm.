function [ X_est,om1,om2 ] = matrixOMP( A1, A2, Y, T, th)
%OMP Summary of this function goes here
%   Detailed explanation goes here
%   Model: Y=A1*X*A2+N
[M1, N1] = size(A1);
[N2, M2] = size(A2);
norm_A1 = diag(A1'*A1);
norm_A2 = diag(A2*A2');
norm_A = norm_A1*norm_A2.';
x_est = zeros(N1*N2,1);
X_est = zeros(N1,N2);

om = [];
om1 = [];
om2 = [];
Aom = [];
for t = 1:T
    R = Y - A1*X_est*A2;
    if (norm(R,'fro')^2<th)
        break;
    end
    MF = (A1'*R*A2')./norm_A;
    [~,idx] = max(abs(MF(1:end)));
    om = [om idx];
    id1 = mod(idx-1,N1) + 1;
    om1 = [om1 id1];
    id2 = floor((idx-1)/N1) + 1;
    om2 = [om2 id2];
    Aom = [Aom, kron(A2(id2,:).',A1(:,id1))];
    x_est(om) = (Aom'*Aom)\(Aom'*Y(:));
    X_est = reshape(x_est, [N1,N2]);
end

end
