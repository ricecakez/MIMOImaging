function [mu,nu,S] = U_ESPRIT2D(Y,M,N,K)
L = size(Y,2);
if mod(M*N,2) == 0
    m = M*N/2;
else
    m = (M*N-1)/2;
    y = Y(m+1,:);
end
Y1 = Y(1:m,:);
Y2 = Y((end-m+1):end,:);
Y21 = IEM(m) * Y2;
if mod(M*N,2) == 0
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        imag(Y1-Y21) real(Y1-Y21)];
else
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        sqrt(2)*real(y) -sqrt(2)*imag(y);
        imag(Y1-Y21) real(Y1-Y21)];
end
[U,D,V] = svd(Z);
Es = U(:,1:K);
J2 = kron(eye(N),[zeros(M-1,1), eye(M-1)]);
J4 = [ zeros((N-1)*M,M),eye((N-1)*M)];
tmp1 = UniMat((M-1)*N)'*J2*UniMat(M*N);
K1 = real(tmp1);
K2 = imag(tmp1);
tmp2 = UniMat((N-1)*M)'*J4*UniMat(M*N);
K3 = real(tmp2);
K4 = imag(tmp2);

Lambda_mu = pinv(K1*Es)*(K2*Es);
Lambda_nu = pinv(K3*Es)*(K4*Es);
[V,D] = eig(Lambda_mu+1i*Lambda_nu);
Xi_mu = real(diag(D));
Xi_nu = imag(diag(D));
mu = atan(Xi_mu)/pi;
nu = atan(Xi_nu)/pi;
A = kr(exp(1i*2*pi*(0:N-1).'*nu.'),exp(1i*2*pi*(0:M-1).'*mu.'));
S = sqrt(abs(diag((A'*A)\A'*(Y*Y'/L)*A/(A'*A))));