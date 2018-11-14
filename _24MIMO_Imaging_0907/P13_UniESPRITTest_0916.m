clear;clc;close all;

SNR = 5;
M = 8;
N = 6;
P = 3;
theta = [10 -8 0]*pi/180;
phi = [20 30 45]*pi/180;
u = sin(theta);
v = sin(phi);
at = exp(1i*pi*(0:M-1).'*u);
ar = exp(1i*pi*(0:N-1).'*v);
L = 50;
S = (randn(P,L) + 1i*randn(P,L))/sqrt(2);
Y = awgn(kr(ar,at)*S,SNR,'measured');
Z = [Y fliplr(eye(M*N))*conj(Y)*fliplr(eye(L))];
QMN = UniMat(M*N);
Q2L = UniMat(2*L);
Gamma = real(QMN'*Z*Q2L);
R = Gamma*Gamma'/2/L;
[V,D] = eig(R);
[Di,ind] = sort(diag(D),'descend');
Vi = V(:,ind);
Es = Vi(:,1:P);

