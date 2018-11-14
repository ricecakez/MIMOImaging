clear;clc;close all;

M = 8;
N = 6;
K = 3;
theta = [10 20 30]*pi/180;
ar = exp(-1i*pi*(0:N-1)'*sin(theta));
at = exp(-1i*pi*(0:M-1)'*sin(theta));
L = 50;
s = real(exp(1i*pi*[0.1;0.2;0.3]*(0:L-1)));
SNR = -5:5:25;
I = 3000;
for nn = 1:length(SNR)
    for i = 1:I
        x = awgn(kr(ar,at)*s,SNR(nn),'measured');
        t_est(i,:) =  sort(RC_ESPRIT(x,M,N,K)*180/pi);
        t_est1(i,:) =  sort(Unitary_ESPRIT_Re(x,M,N,K)*180/pi);
    end
    RMSE(nn) = sqrt(sum(sum((t_est - theta).^2)))/K/I;
end
plot(SNR,RMSE)
% figure
% plot(1:I,t_est,'.','linewidth',100)