function theta_est = Unitary_ESPRIT_Re(x,M,N,K)
G0 = [eye(M) zeros(M,N-1)];
G = [];
for m = 0:N-1
    G = [G;circshift(G0,m,2)];
end
W = G'*G;
%½µÎ¬±ä»»
Y = W^(-0.5)*G'*x;
[P,L] = size(Y);
% Y2 = Y(2:end,:);
Y3 = conj(Y);
Z = [IEM(P)*Y;Y3];
Q = UniMat(2*P);
Zc = real(Q'*Z);
Rc = Zc*Zc'/L;
Wz = diag([diag(W).^(-0.5);diag(W).^(-0.5)]);
J1 = [eye(P) zeros(P)];
J2 = [zeros(P) eye(P)];
K1 = Re(Q'*J2*Wz*Q);