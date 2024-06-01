function hn = channel_norm(h)
[N, K] = size(h);

hn = zeros(size(h));

for k = 1:K
    hk = h(:, k);
    hn(:, k) = hk / sqrt(sum(abs(hk(:)).^2)) * sqrt(N );
end


end

