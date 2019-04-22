function [A,B,C] = TALS(X,R)
Xnorm = norm(X);
sz = size(X);
% A0 = randn(sz(1),R);
% B0 = randn(sz(2),R);
% C0 = randn(sz(3),R);
% A0 = exp(1i*2*pi*(0:sz(1)-1).'*rand(1,R));
% B0 = exp(1i*2*pi*(0:sz(2)-1).'*rand(1,R));
% C0 = exp(1i*2*pi*(0:sz(3)-1).'*rand(1,R));
% lambda = randn(R,1);
X1 = double(tenmat(X,1));
[A0,~,~] = svds(X1,R);
X2 = double(tenmat(X,2));
[B0,~,~] = svds(X2,R);
X3 = double(tenmat(X,3));
[C0,~,~] = svds(X3,R);

% U{1,1} = A0;
% U{1,2} = B0;
% U{1,3} = C0;
% lambda = ones(3,1);
Y = full(ktensor({A0,B0,C0}));
Res0 = 0;
% Res = sqrt(Xnorm^2+norm(Y)^2-2 * innerprod(X,Y));
Res1 = norm(X-Y);
dRes = abs(Res1-Res0);
ResTh = 1e-4;
Niter = 200;
ii = 1;
A = A0;
B = B0;
C = C0;
X1 = double(tenmat(X,1));
X2 = double(tenmat(X,2));
X3 = double(tenmat(X,3));
while (dRes(ii) >= ResTh) && (ii<=Niter)
    C = X3/(krb(B,A)).';
    B = X2/(krb(C,A)).';    
    A = X1/(krb(C,B)).';
%     lambda = ones(R,1);
%     for i= 1:R
%         lambda(i) = A(1,i)*B(1,i)*C(1,i);
%         A(:,i) = A(:,i)/A(1,i);
%         B(:,i) = B(:,i)/B(1,i);
%         C(:,i) = C(:,i)/C(1,i);
%     end
    Y = full(ktensor({A,B,C}));
    Res0 = Res1;
%     Res = abs(sqrt(Xnorm^2+norm(Y)^2-2 * innerprod(X,Y)));
    Res1 = norm(X-Y);
    ii = ii + 1;
    dRes(ii) = abs(Res1-Res0);
end
end