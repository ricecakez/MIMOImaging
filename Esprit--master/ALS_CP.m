function [A,lambda] = ALS_CP(X,R)
N = ndims(X);
Xnorm = norm(X);
lambda = ones(R,1);
for n = 1:N
    x = double(tenmat(X,n));
    [U,D,~] = svd(x);
    A{n} = U(:,1:R);
    for r = 1:R
        lambda(r) = lambda(r)*sqrt(D(r,r));
        %         A{n}(:,r) = A{n}(:,r)/norm(A{n}(:,r));
    end
end
U = ktensor(lambda,A);
Unorm = norm(U);
Res0 = 0;
Normresidual = sqrt(Xnorm^2 + Unorm^2 - 2*innerprod(X,U));
Niter = 100;
DresTh = 1e-4;
ii = 1;
dres = Normresidual - Res0;
for n = 1:N
    if ~isempty(A{n})
        AtA(:,:,n) = A{n}'*A{n};
    end
end
while (ii <= Niter) && (dres >DresTh)
    Res0 = Normresidual;
    for n = 1:N
        V = prod(AtA(:,:,[1:n-1 n+1:N]),3);
        Anew = mttkrp(X,A,n);
        Anew = Anew / V;
        lambda = max( max(abs(Anew),[],1), 1 )'; %max-norm
        Anew = bsxfun(@rdivide, Anew, lambda');
        A{n} = Anew;
        AtA(:,:,n) = A{n}'*A{n};
    end
    U = ktensor(lambda,A);
    Unorm = norm(U);
    Normresidual = sqrt(Xnorm^2 + Unorm^2 - 2*innerprod(X,U));
    dres = abs(Normresidual - Res0);
    ii = ii+1;
end