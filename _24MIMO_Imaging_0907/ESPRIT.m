function f = ESPRIT(x,M,K)
%% spatial smoothing & estimation of the correlation matrix
% x1 = x*fliplr(eye(length(x)));
% xs = buffer(x,M,M-1,'nodelay');
L = length(x);
rx = xcorr(x(1:end),M,'unbiased');
% rxy = xcorr(x(1:end-1),x(2:end),'biased');
for m = 1:M
    Rxx(m,:) = rx((M+m):-1:(M+m-M+1));
    Rxy(m,:) = rx((M+m-1):-1:(M+m-M));
end
% xs = flipud(xs);
% L = size(xs,2);
% % Rxx = xs(:,1:end-1)*xs(:,1:end-1)'/(L-1);
% Rxy = xs(:,1:end-1)*xs(:,2:end)'/(L-1);
%% SVD and Calculate the minimum eigenvalue
[V,D] = svd(Rxx);
ev = diag(D); emin = ev(end);
Z = zeros(M);
for m = 2:M
    Z(m,m-1) = 1;
end
Cxx = Rxx - emin*eye(M);
Cxy = Rxy - emin*Z;
% [V,D] = eig(Cxx,Cxy);
[U,D,V] = svd(Cxx);
D1 = D(1:K,1:K);
U1 = U(:,1:K);
V1 = V(:,1:K);
Cxy1 = U1'*Cxy*V1;
[V,D] = eig(D1,Cxy1);
z = diag(D);
f = angle(z)/2/pi;


