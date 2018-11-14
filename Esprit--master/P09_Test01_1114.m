clear;clc;close all;
%测试角度兼并问题
t = [0.1 0.1 0.3];
r = [0.2 0.3 0.4];
M = 8;
N = 8;
K = length(t);
L = 200;
At = exp(1i*2*pi*(0:M-1).'*t);
Ar = exp(1i*2*pi*(0:M-1).'*r);
A = kr(Ar,At);
A1 = A(1:K,:);
dd = rank(A1);
% Phir = diag(exp(1i*2*pi*r));
% AA = [];
% for n = 1:N
%     AA = [AA;At*Phir^(n-1)];
% end
% rank(AA)
B = exp(1i*2*pi*randn(K,L));
X = A*B;
[t1,r1] = PM2D(X,M,N,K);
[t1,r1] = sort2D(t1,r1,1);
[t1,r1] = ESPRIT_2D_Modi(X,M,N,K);
[t1,r1] = sort2D(t1,r1,1)
