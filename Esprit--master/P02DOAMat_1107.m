clear;clc;close all;

K = 3;
alpha = [60,70,80]*pi/180;
beta = [60,76,80]*pi/180;
L = 100;
M = 8;
N = 8;
Ar = exp(1i*pi*(0:N-1).'*cos(beta));
At = exp(1i*pi*(0:M-1).'*cos(alpha));
A = kr(Ar,At);
d = exp(1i*2*pi*rand(K,1)*(0:L-1));
P = d*d'/L;
% A'*A
% P'*P;
SNR = 0;
x0 = A*d;
P1 = x0*x0'/L;
P2 = A*P*A';
x = x0;
% x = awgn(x0,SNR,'measured');

x1 = x(1:M*(N-1),:);
x2 = x(end-M*(N-1)+1:end,:);
R1 = x1*x1'/L;
R2 = x2*x2'/L;
[V,D] = eig(R1);
[ds,inds] = sort(diag(D),'descend');
sigma2 = mean(ds(K+1:end))';
% L = diag(ds(1:K));
% U = V(:,inds(1:K));
Sigma = [zeros((N-2)*M,M), eye((N-2)*M)
    zeros(M), zeros(M,(N-2)*M)];
G1 = R1 - sigma2*eye((N-1)*M);
G2 = R2 - sigma2*eye((N-1)*M);
[V,D] = eig(G1);
[ds,inds] = sort(diag(D),'descend');
% sigma2 = mean(ds(K+1:end))';
L = diag(ds(1:K));
U = V(:,inds(1:K));
% G11 = U*L*U';
G11 = zeros((N-1)*M);
for k = 1:K
    G11 = G11 + 1/L(k,k)*U(:,k)*U(:,k).';
end
R = G2*G11;
[A,Phi] = eig(R);
[phi1,inds] = sort(diag(Phi),'descend');
alpha_est = acos(angle(phi1(1:K))/pi)*180/pi;
A1 = A(:,inds(1:K));
