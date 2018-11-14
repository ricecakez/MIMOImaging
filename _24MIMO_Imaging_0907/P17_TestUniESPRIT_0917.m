clear;clc;close all;

P = 2;
N = 50;
M = 16;
S = exp(1i*rand(P,N));
theta = [10 20]*pi/180;
SNR = 5;
a = exp(1i*pi*(0:M-1).'*sin(theta));
x = awgn(a*S,SNR,'measured');
tic
Rxx = x*x'/N;
[V,D] = eig(Rxx);
[di,inds] = sort(diag(D),'descend');
Vi = V(:,inds);
Us = Vi(:,1:P);
U1 = Us(1:end-1,:);
U2 = Us(2:end,:);
Psi = inv(U1'*U1)*U1'*U2;
[V,D] = eig(Psi);
d = diag(D);
theta_est = asin(angle(d)/pi)*180/pi

toc

QM = UniMat(M);
QN = UniMat(2*N);
PIM = IEM(M);
PIN = IEM(N);

tic
z = [x PIM*conj(x)*PIN];
G = real(QM'*z*QN);
[U,D,V]=svd(G);
Us = QM*U(:,1:P);
U1 = Us(1:end-1,:);
U2 = Us(2:end,:);
Psi = inv(U1'*U1)*U1'*U2;
[V,D] = eig(Psi);
d = diag(D);
theta_est = asin(angle(d)/pi)*180/pi
toc



