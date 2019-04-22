function Z = UniTrans(X)
Xc = tensor(conj(double(X)));
dims = size(X);
M = dims(1);
N = dims(2);
P = dims(3);
X1 = ttm(Xc,{IEM(M),IEM(N),[zeros(P-1,1) IEM(P-1)]});
Y = tensor(cat(3,double(X1),double(X)));
P1 = size(Y,3);
Z = tensor(real(double(ttm(Y,{UniMat(M)',UniMat(N)',UniMat(P1)'}))));
end