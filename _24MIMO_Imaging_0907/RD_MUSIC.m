function [t,r] = RD_MUSIC(X,f2,M,N,K)
L = size(X,2);
Rx = X*X'/L;

e1 = [1;zeros(M-1,1)];
[V,D] = eig(Rx);
[d1,inds] = sort(diag(D),'descend');
En = V(:,inds(K+1:end));
for i = 1:length(f2)
    ar = exp(1i*2*pi*(0:N-1).'*f2(i));
    Q = (kron(ar,eye(M)))'*En*En'*(kron(ar,eye(M)));
    Pr(i) = e1'/Q*e1;
end

indr = find_peak(Pr,K);
r = f2(indr);
P = [ones(M,1) [0:M-1].'];
for k = 1:K
    ar = exp(1i*2*pi*(0:N-1).'*r(k));
    Q = (kron(ar,eye(M)))'*En*En'*(kron(ar,eye(M)));
    at = (Q\e1)/Pr(indr(k));
    gt = phase(at)/2/pi;
    C = (P.'*P)\P.'*gt;
    t(k) = C(2);
end