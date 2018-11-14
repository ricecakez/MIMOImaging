clear;clc;close all;

K = 3;
alpha = [60,70,80]*pi/180;
beta = [60,76,80]*pi/180;
L = 100;
M = 8;
N = 8;
Ar = exp(1i*pi*(0:N-1).'*sin(beta));
At = exp(1i*pi*(0:M-1).'*sin(alpha));
A = kr(Ar,At);
B = exp(1i*2*pi*rand(K,L));
SNR = 0;
X = awgn(A*B,SNR,'measured');
dt = 0.1;
alpha_est = (-90:dt:90-dt);
beta_est = (-90:dt:90-dt);

Rx = X*X'/L;
for p = 1:length(alpha_est)
    for q = 1:length(beta_est)
        ar = exp(1i*pi*(0:N-1).'*sin(beta_est(q)*pi/180));
        at = exp(1i*pi*(0:N-1).'*sin(alpha_est(p)*pi/180));
        f(p,q) = 1/((kron(ar,at))'*inv(Rx)*(kron(ar,at)));
    end
end
figure
imagesc(beta_est,alpha_est,abs(f))