clear;clc;close all;

% t = [10 20 30];
% r = [20 30 40];
t = [0.1 0.2 0.3];
r = [0.2 0.3 0.4];
K = 3;
M = 8;
N = 8;
L = 200;
At = exp(1i*2*pi*(0:M-1).'*t);
Ar = exp(1i*2*pi*(0:N-1).'*r);
A = kr(Ar,At);
B = exp(1i*randn(K,L));
S = mean(abs(B).^2,2);
X0 = A*B;
SNR = -5:5:25;
MM = 300;
Re = zeros(size(SNR));
R1 = zeros(size(SNR));
R11 = zeros(size(SNR));
R2 = zeros(size(SNR));
R3 = zeros(size(SNR));
f2 = 0:0.01:0.5;
for nn = 1:length(SNR)
for mm = 1:MM
X = awgn(X0,SNR(nn),'measured');
tic
[te,re] = RD_MUSIC(X,f2,M,N,K);
toc
[te,re] = sort2D(te,re,1);
Re(nn) = Re(nn) + mean((te - t).^2 + (re - r).^2);
tic
[t1,r1] = PM2D(X,M,N,K);
toc
[t1,r1] = sort2D(t1,r1,1);
R1(nn) = R1(nn) + mean((t1.' - t).^2 + (r1.' - r).^2);

tic
[t11,r11] = PM2D1(X,M,N,K);
toc
[t11,r11] = sort2D(t11,r11,1);
R11(nn) = R11(nn) + mean((t11.' - t).^2 + (r11.' - r).^2);

tic
[t2,r2,s2] = ESPRIT_2D(X,M,N,K);
toc
[t2,r2] = sort2D(t2,r2,1);
R2(nn) = R2(nn) + mean((t2.' - t).^2 + (r2.' - r).^2);

tic
[t3,r3,s3] = Unitary_ESPRIT_2D1(X,M,N,K);
toc
[t3,r3] = sort2D(t3,r3,1);
R3(nn) = R3(nn) + mean((t3.' - t).^2 + (r3.' - r).^2);
% ta1 = -asin(t1)/pi*180;
% ra1 = -asin(r1)/pi*180;
% figure(nn)
% scatter(t1,r1,'.','k');
% hold on
% scatter(te,re,'.','p');
% scatter(t1,r1,'.','g');
% scatter(t2,r2,'.','r');
% scatter(t3,r3,'.','b');
end
Re(nn) = sqrt(Re(nn)/M);
R1(nn) = sqrt(R1(nn)/M);
R11(nn) = sqrt(R11(nn)/M);
R2(nn) = sqrt(R2(nn)/M);
R3(nn) = sqrt(R3(nn)/M);
end
figure
plot(SNR,mag2db(Re),SNR,mag2db(R1),SNR,mag2db(R11),SNR,mag2db(R2),SNR,mag2db(R3))