function U = CP_ALS(X,K,iterNum,thresh)
N = ndims(X);
normX = norm(X);
for n = 1:N
    Xn = double(tenmat(X,n));
    [U,~,~] = svds(Xn,K);
    U_ini{n} = U(:,1:K);
end
Y = ktensor(U_ini);
res0 = sqrt(normX^2 + norm(Y)^2 - 2 * innerprod(X,Y));
dres = abs(res0);
ii = 1;
U = U_ini;
while (dres(ii) >= thresh) && (ii<=iterNum)
    for n = 1:N
        Xn = double(tenmat(X,n));
        Vn = U_kr(U,K,n);
        U{n} = Xn/Vn.';
    end
    Y = ktensor(U);
    res = sqrt(normX^2 + norm(Y)^2 - 2 * innerprod(X,Y));
    ii = ii + 1;
    dres(ii) = abs(res - res0);
    res0 = res;
end