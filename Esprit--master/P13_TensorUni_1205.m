clear;clc;close all;

M = 9;
K = 3;
N = 9;
phi = [10 20 30 40 50 60];
theta = [20 30 40 50 60 70];
fd = -0.2:0.1:0.2;
% hold on
alpha = -sin(phi*pi/180)/2;
beta = -sin(theta*pi/180)/2;
P = 10;
At = exp(1i*2*pi*(-(M-1)/2:(M-1)/2).'*alpha(1:K));
Ar = exp(1i*2*pi*(-(N-1)/2:(N-1)/2).'*beta(1:K));
S = exp(1i*2*pi*(0:P-1).'*fd(1:K));
SNR = 5;
X0 = ttm(tensorI(K),{At,Ar,S});
X1 = double(tenmat(X0,3)).';
X2 = awgn(X1,SNR,'measured');
X = tensor(reshape(X2,[M,N,P]));

Z0 = UniTrans(X);
% Xc = tensor(conj(double(X0)));
% X1 = ttm(Xc,{IEM(M),IEM(N),[zeros(P-1,1) IEM(P-1)]});
% % X2 = ttm(tensorI(K),{IEM(M)*conj(At),IEM(N)*conj(Ar),IEM(P)*conj(S)});
% Y = tensor(cat(3,double(X1),double(X0)));
% Z = ttm(Y,{UniMat(M)',UniMat(N)',UniMat(2*P-1)'});
% S1 = real(UniMat(2*P)'*(cat(1,S,IEM(P)*conj(S))));
% At_r = real(UniMat(M)'*IEM(M)*conj(At));
% Jt1 = [eye(M-1) zeros(M-1,1)];
Jt2 = [zeros(M-1,1) eye(M-1)];
tmp = UniMat(M-1)'*Jt2*UniMat(M);
Kt1 = real(tmp);
Kt2 = imag(tmp);
Jr2 = [zeros(N-1,1) eye(N-1)];
tmp = UniMat(N-1)'*Jr2*UniMat(N);
Kr1 = real(tmp);
Kr2 = imag(tmp);

Jd2 = [zeros(2*P-2,1) eye(2*P-2)];
tmp = UniMat(2*P-2)'*Jd2*UniMat(2*P-1);
Kd1 = real(tmp);
Kd2 = imag(tmp);

% phi_t = pinv(Kt1*At_r)*Kt2*At_r;

% Z1 = tensor(real(double(Z)));
[A,B,C] = TALS(Z0,K);
phi_a = pinv(Kt1*A)*Kt2*A;
asin(-atan(diag(phi_a))/pi*2)*180/pi
phi_b = pinv(Kr1*B)*Kr2*B;
asin(-atan(diag(phi_b))/pi*2)*180/pi
phi_d = pinv(Kd1*C)*Kd2*C;
atan(diag(phi_d))/pi
