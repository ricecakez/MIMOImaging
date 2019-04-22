function V = U_kr(U,R,n)
N = length(U);
nn = [1:n-1 n+1:N];
V = ones(1,R);
for ii = 1:length(nn)
    V = kr(U{nn(ii)},V);
end