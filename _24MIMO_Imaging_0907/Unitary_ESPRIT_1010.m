function [f1,f2,S] = Unitary_ESPRIT_1010(x,M,N,d)
N0 = N-1;
X = buffer(x,M*N0,M*N0-1,'nodelay');
% Y = [X IEM(M)*conj(x1)*IEM(N)];
% if mod(M,2) == 0
    m = M*N0/2;
    X1 = X(1:m,:);
    X2 = X(m+1:end,:);
    T1 = real(X1+IEM(m)*X2);
    T2 = -imag(X1+IEM(m)*X2);
    T3 = imag(X1-IEM(m)*X2);
    T4 = real(X1-IEM(m)*X2);
    Z = [T1 T2
        T3 T4];
% else
%     m = (M-1)/2;
%     X1 = X(1:m,:);
%     X2 = X(m+1,:);
%     X3 = X(m+2:end,:);
%     T1 = real(X1+IEM(m)*X3);
%     T2 = -imag(X1+IEM(m)*X3);
%     T3 = sqrt(2)*real(X2);
%     T4 = -sqrt(2)*imag(X2);
%     T5 = imag(X1-IEM(m)*X3);
%     T6 = real(X1-IEM(m)*X3);
%      Z = [T1 T2
%          T3 T4
%         T5 T6];
% end
[U,D,V] = svd(Z);
Es = U(:,1:d);
J2 = (UniMat(N0*(M-1)))'*kron(eye(N0),[zeros(M-1,1) eye(M-1)])*UniMat(N0*M);
K1 = real(J2);
K2 = imag(J2);
Psi1 = pinv(K1*Es)*K2*Es;
J4 = (UniMat(M*(N0-1)))'*kron([zeros(N0-1,1) eye(N0-1)],eye(M))*UniMat(N0*M);
K3 = real(J4);
K4 = imag(J4);
Psi2 = pinv(K3*Es)*K4*Es;
[T,Phi] = eig(Psi1+1i*Psi2); 
f1 = atan(real(diag(Phi)))/pi;
f2 = atan(imag(diag(Phi)))/pi;
S = abs(T\Es'*(UniMat(M*N0))'*x(1:M*N0));
end
