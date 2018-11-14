clear;clc;close all;

load('RxSignal_plane_1106.mat');
SNR = 5;
[phi,psi] = sort2D(phi,psi,1);
A1 = exp(1i*2*pi*(0:N0-1).'*phi.');
A2 = exp(1i*2*pi*(0:Nv-1).'*psi.');
A = kr(A2,A1);
S = 12 * exp(1i*2*pi*rand(I,1)*(0:K-1));
X0 =  A*S;
X = awgn(X0,SNR,'measured');
tic
[p1,s1] = Unitary_ESPRIT_2D1(X,N0,Nv,I);
toc
[p1,s1] = sort2D(p1,s1,1);
tic
[p2,s2] = PM2D(X,N0,Nv,I);
toc
[p2,s2] = sort2D(p2,s2,1);
f2 = -0.2:0.001:0.2;
tic
[p3,s3] = RD_MUSIC(X,f2,N0,Nv,I);
toc
[p3,s3] = sort2D(p3.',s3.',1);
