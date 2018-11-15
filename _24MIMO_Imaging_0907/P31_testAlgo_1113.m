clear;clc;close all;

load('data\RxSignal_plane_1106.mat');
SNR = 5;
[phi,psi] = sort2D(phi,psi,1);
u = phi*N0*d0;
v = psi*lambda/2*R0/d;
figure
scatter(u,v,12)
A1 = exp(1i*2*pi*(0:N0-1).'*phi.');
A2 = exp(1i*2*pi*(0:Nv-1).'*psi.');
A = kr(A2,A1);
% A*A'
S = 12 * exp(1i*2*pi*rand(I,1)*(0:K-1));
P = S*S'/K;
% diag(P)
X0 =  A*S;
% X = X0;
X = awgn(X0,SNR,'measured');
tic
[p1,s1,S1] = Unitary_ESPRIT_2D(X,N0,Nv,I);
toc
tic
[p2,s2] = Unitary_ESPRIT_2D1115(X,N0,Nv,I);
toc

ms = sort_matrix([p1 s1 S1],'ascend',1);
p1 = ms(:,1);
s1 = ms(:,2);
S1 = ms(:,3);

ms = sort_matrix([p2 s2],'ascend',1);
p2 = ms(:,1);
s2 = ms(:,2);
% S2 = ms(:,3);

u_est = p1*N0*d0;
v_est = s1*lambda/2*R0/d;
figure
scatter(u_est,v_est,S1);
x_est = v_est - u_est*u0(1);
y_est =  (- x_est*u0(1)-u_est)/u0(2);
figure
scatter(x_est,y_est,S1);

u_est2 = p2*N0*d0;
v_est2 = s2*lambda/2*R0/d;
figure
scatter(u_est2,v_est2,12);
x_est2 = v_est2 - u_est2*u0(1);
y_est2 =  (- x_est2*u0(1)-u_est2)/u0(2);
figure
scatter(x_est2,y_est2,S1);
% tic
% [p2,s2] = PM2D(X,N0,Nv,I);
% toc
% [p2,s2] = sort2D(p2,s2,1);
% f2 = -0.2:0.001:0.2;
% tic
% [p3,s3] = RD_MUSIC(X,f2,N0,Nv,I);
% toc
% [p3,s3] = sort2D(p3.',s3.',1);
