clear;clc;close all;

K = 3;
alpha = [10,20,30]*pi/180;
beta = [20,30,40]*pi/180;
L = 100;
M = 8;
N = 8;
Ar = exp(-1i*pi*(0:N-1).'*sin(beta));
At = exp(-1i*pi*(0:M-1).'*sin(alpha));
A = kr(Ar,At);
B = randn(K,L);
% B = exp(1i*2*pi*rand(K,L));
SNR = 10;
X0 = A*B;
P = 40;
for p = 1:P
X = awgn(X0,SNR,'measured');
Rx = X*X'/L;
R1 = inv(Rx);
dt = 0.1;
beta1 = -90:dt:(90-dt);
e1 = [1;zeros(M-1,1)];
for l = 1:length(beta1)
    ar = exp(-1i*pi*(0:N-1).'*sin(beta1(l)*pi/180));
    Q = (kron(ar,eye(M)))'*R1*(kron(ar,eye(M)));
    f(l) = e1'/Q*e1;
end
f = abs(f);
% figure
% plot(beta1,f);
beta_est = beta1(find_peak(f,K));
P1 = [ones(M,1),pi*(0:M-1).'];
for k = 1:K
    ar = exp(-1i*pi*(0:N-1).'*sin(beta_est(k)*pi/180));
    Q = (kron(ar,eye(M)))'*R1*(kron(ar,eye(M)));
    at = Q\e1/(e1'/Q*e1);
    gk = -phase(at);
    ck = inv(P1.'*P1)*P1.'*gk;
    alpha_est(k)=asin(ck(2))/pi*180;
end
figure(1)
scatter(alpha_est,beta_est,'.','k')
hold on
end