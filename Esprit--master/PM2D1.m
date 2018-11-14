function [t1,r1] = PM2D1(X,M,N,K)
L = size(X,2);
Rx = X*X'/L;

G = Rx(:,1:K);
H = Rx(:,K+1:end);
Pc = ((G'*G)\G'*H)';
P = [eye(K);Pc];
Q = P*(P'*P)^(-1/2);
P1 = Q(1:end-M,:);
P2 = Q(M+1:end,:);
Psir = (P1'*P1)\P1'*P2;
[Vr,Dr] = eig(Psir);
r1 = angle(diag(Dr))/2/pi;

Te = rowPermutateMat(M,N);
Pr = Te*Q;
B = Pr*Vr;
B1 = B(1:end-N,:);
B2 = B(N+1:end,:);
Phit = B1\B2;
t1 = angle(diag(Phit))/2/pi;