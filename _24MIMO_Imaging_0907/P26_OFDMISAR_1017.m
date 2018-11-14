clear;clc;close all;

load('OFDMISAREcho_mat');

SNR = 5;
Y0 = awgn(X,SNR,'measured');
for k = 1:K
    t = T_min + (k-1)*T + tc + (0:1/fs:(tb-1/fs));
    t = t(:);
    Y(:,k) = Y0((L_c+1):L_T,k).*exp(-1i*2*pi*f0*t);
    Y1(:,k) = fft(Y(:,k))./a(:,k)/N;
end 
