function RS = ssp(R,K)
[M,N] = size(R);
L = M-K+1;
RS = zeros(K);
for l = 1:L
    RS = RS + R(l:l+K-1,l:l+K-1);
end
RS = RS/N;