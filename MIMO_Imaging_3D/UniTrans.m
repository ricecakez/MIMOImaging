function Z = UniTrans(X,dim)
Xc = tensor(conj(double(X)));
dims = size(X);
M = dims(1);
N = dims(2);
P = dims(3);
if dim == 1
    X1 = ttm(Xc,{[zeros(M-1,1) IEM(M-1)],IEM(N),IEM(P)});
else if dim == 2
        X1 = ttm(Xc,{IEM(M),[zeros(N-1,1) IEM(N-1)],IEM(P)});
    else
        X1 = ttm(Xc,{IEM(M),IEM(N),[zeros(P-1,1) IEM(P-1)]});
    end
end
Y = tensor(cat(dim,double(X1),double(X)));
M1 = size(Y,1);
N1 = size(Y,2);
P1 = size(Y,3);
Z = tensor(real(double(ttm(Y,{UniMat(M1)',UniMat(N1)',UniMat(P1)'}))));
end