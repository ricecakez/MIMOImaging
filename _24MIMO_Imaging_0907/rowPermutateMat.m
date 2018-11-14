function Te = rowPermutateMat(M,N)
I = eye(M*N);
p = [];
for m = 1:M
    p = [p m:M:M*N];
end
Te = I(p,:);