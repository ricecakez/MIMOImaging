clear;clc;close all;

K = 3;
theta = [10,20,30,40,50,60,70];
phi = [20,15,40,50,60,70,80];
figure(1)
scatter(phi(1:K),theta(1:K),'ro')
hold on
L = 50;
M = 9;
N = 9;
rcs = rand(L,K);
fd = 0.1*(-(K-1)/2:(K-1)/2);
At = exp(-1i*pi*(0:M-1).'*sin(theta(1:K)*pi/180));
Ar = exp(-1i*pi*(0:N-1).'*sin(phi(1:K)*pi/180));
B = exp(1i*2*pi*(0:L-1).'*fd);
X0 = kr(Ar,At)*B.';
SNR = 10;
MM = 500;
for ii = 1:MM
X = awgn(X0,SNR,'measured');
X1 = reshape(X,[M,N,L]);
% p = cp_als(tensor(X1),K,'dimorder',[3 2 1]);
% [At1,Ar1,B1] = TALS(tensor(X1),K);
% U = CP_ALS(tensor(X1),K,100,1e-4);
% At1 = double(U{1});
% Ar1 = double(U{2});
% B1 = double(U{3});
% % for n = 1:N
% %     Y(:,n) = reshape(X((n-1)*M+(1:M),:).',[M*L,1]);
% % end
% % for m = 1:M
% %     Z(:,m) = reshape(X(m:M:end,:),[N*L,1]);
% % end
% % At1 = exp(-1i*pi*(0:M-1).'*rand(1,K));
% % Ar1 = exp(-1i*pi*(0:N-1).'*rand(1,K));
% % B1 = (pinv(kr(Ar1,At1))*X).';
% % Ar1 = (pinv(kr(At1,B1))*Y).';
% % At1 = (pinv(kr(B1,Ar1))*Z).';
% % SSR(1) = sqrt(sum(sum((abs(X - (kr(Ar1,At1)*B.'))).^2)));
% % DSSR(1) = 20*log10(SSR);
% % i = 1;
% % while (DSSR(i) > -100)
% %     i = i + 1;
% %     B1 = (pinv(kr(Ar1,At1))*X).';
% %     Ar1 = (pinv(kr(At1,B1))*Y).';
% %     At1 = (pinv(kr(B1,Ar1))*Z).';
% %     SSR(i) = sqrt(sum(sum((abs(X - (kr(Ar1,At1)*B.'))).^2)));
% %     DSSR(i) = 20*log10(abs(SSR(i)-SSR(i-1)));
% %     if abs(DSSR(i)-DSSR(i-1))<0.001
% %         break;
% %     end
% % end
% P1 = [ones(N,1) (0:(N-1)).'];
% P2 = [ones(M,1) (0:(M-1)).'];
% P3 = [ones(L,1) (0:(L-1)).'];
% for k = 1:K
%     h = -phase(Ar1(:,k))/pi;
%     g = -phase(At1(:,k))/pi;
%     d = phase(B1(:,k))/2/pi;
%     c1 = (P1.'*P1)\P1.'*h;
%     phi_est(k) = asin(c1(2,:))*180/pi;
%     c2 = (P2.'*P2)\P2.'*g;
%     theta_est(k) = asin(c2(2,:))*180/pi;
%     c3 = (P3.'*P3)\P3.'*d;
%     fd_est(k) = c3(2);
% end
f = spect_est(tensor(X1),K,100,1e-5);
theta_est = asin(-f(1,:)*2)*180/pi;
phi_est = asin(-f(2,:)*2)*180/pi;
fd_est = f(3,:);
figure(1)
scatter(phi_est,theta_est,'k+')
grid on
hold on
xlim([-90,90])
ylim([-90,90])
clear SSR;clear DSSR;
end