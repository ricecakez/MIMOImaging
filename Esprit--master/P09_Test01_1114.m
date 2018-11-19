clear;clc;close all;
%测试角度兼并问题
t = [0.1 0.2 0.3];
r = [0.2 0.3 0.3];
M = 8;
N = 8;
K = length(t);
L = 200;
At = exp(1i*2*pi*(0:M-1).'*t);
Ar = exp(1i*2*pi*(0:M-1).'*r);
A = kr(Ar,At);
eye(M*N)-A/(A'*A)*A';
% for k = 1:K
%     d1 = UniMat(N)'*Ar(:,k)*At(:,k).'*conj(UniMat(M));
%     tmp = UniMat(M*N)'*A(:,k);
%     d2 = reshape(tmp,[N,M]);
% end
    

% A1 = A(1:K,:);
% dd = rank(A1);
% Phir = diag(exp(1i*2*pi*r));
% AA = [];
% for n = 1:N
%     AA = [AA;At*Phir^(n-1)];
% end
% rank(AA)
B = 12*exp(1i*2*pi*randn(K,L));
X = A*B;
SNR = 0;
X = awgn(X,SNR,'measured'); 
% [t1,r1] = PM2D(X,M,N,K);
% [t1,r1] = sort2D(t1,r1,1);
tic
[t2,r2,s2] = U_ESPRIT2D(X,M,N,K);
toc
figure
scatter(t2,r2,s2)
% [t2,r2] = sort2D(t2,r2,1)
