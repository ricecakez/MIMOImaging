function theta_est = RC_ESPRIT(x,M,N,K)
G0 = [eye(M) zeros(M,N-1)];
G = [];
for m = 0:N-1
    G = [G;circshift(G0,m,2)];
end
W = G'*G;
%½µÎ¬±ä»»
y = W^(-0.5)*G'*x;
L = size(x,2);
Ry = y*y'/L;
[V,D] = eig(Ry);
[ds,ind] = sort(diag(D),'descend');
Vs = V(:,ind);
E0 = Vs(:,1:K);
E = W^(-0.5)*E0;
E1 = E(1:M+N-2,:);
E2 = E(2:end,:);
Psi = pinv(E1)*E2;
[V1,D1] = eig(Psi);
theta_est = asin(-angle(diag(D1))/pi);