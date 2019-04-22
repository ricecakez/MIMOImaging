function I = tensorI(N)
I0 = zeros(N,N,N);
for n = 1:N
    I0(n,n,n) = 1;
end
I = tensor(I0);