clear;clc;close all;

M = 4;
K = 2;
N = 6;
phi = [10 20 30 40 50 60];
theta = [20 30 40 50 60 70];
fd = -0.3:0.1:0.3;
figure
scatter(phi(1:K),theta(1:K),'*')
hold on
alpha = -sin(phi*pi/180)/2;
beta = -sin(theta*pi/180)/2;
P = 20;
At = exp(1i*2*pi*(-(M-1)/2:(M-1)/2).'*alpha(1:K));
% IEM(M)*conj(At)-At
Ar = exp(1i*2*pi*(-(N-1)/2:(N-1)/2).'*beta(1:K));
% IEM(N)*conj(Ar)-Ar
S = exp(1i*2*pi*(0:P-1).'*fd(1:K));
SNR = [-20:5:20];
rmse = zeros(size(SNR));
rmse1 = zeros(size(SNR));
rmse2 = zeros(size(SNR));
X0 = kr(Ar,At)*S.';
MM = 500;
J11 = [eye(M-1) zeros(M-1,1)];
J12 = [zeros(M-1,1) eye(M-1)];
J21 = [eye(N-1) zeros(N-1,1)];
J22 = [zeros(N-1,1) eye(N-1)];
Jm1 = kron(eye(N),J11);
Jm2 = kron(eye(N),J12);
Jn1 = kron(J21,eye(M));
Jn2 = kron(J22,eye(M));
P1 = [ones(N,1) (0:(N-1)).'];
P2 = [ones(M,1) (0:(M-1)).'];
P3 = [ones(P,1) (0:(P-1)).'];
for nn = 1:length(SNR)
for mm = 1:MM
X1 = awgn(X0,SNR(nn),'measured');
% R = X0*X0';
% [V,D] = eig(R);
% [di,inds] = sort(diag(D),'descend');
% Es = V(:,inds(1:K));
[f1,f2] = ESPRIT_2D(X1,M,N,K);
phi1 = asin(-f1*2)*180/pi;
theta1 = asin(-f2*2)*180/pi;
% figure(1)
% scatter(phi1,theta1,'s')
X2 = reshape(X1,[M,N,P]);
X = tensor(X2);
% X = ttm(tensorI(K),{At,Ar,S});
% T = Uni_HOSVD(X);
T = hosvd(X,1,'ranks',[min(M,K),min(N,K),K],'sequential',false);
tic
[f1,f2,f3] = uni_spect_est(X,K);
toc
phi_est1 = asin(-f1*2).'*180/pi;
theta_est1 = asin(-f2*2).'*180/pi;
fd_est = f3;
% tic
% [At1,Ar1,B1] = TALS(X,K);
% toc
% % At1 = CP.U{1};
% % Ar1 = CP.U{2};
% % B1 = CP.U{3};
% for k = 1:K
%     h = -phase(Ar1(:,k))/pi;
%     g = -phase(At1(:,k))/pi;
%     d = phase(B1(:,k))/2/pi;
%     c1 = (P1.'*P1)\P1.'*h;
%     theta_est1(k,1) = asin(c1(2,:))*180/pi;
%     c2 = (P2.'*P2)\P2.'*g;
%     phi_est1(k,1) = asin(c2(2,:))*180/pi;
%     c3 = (P3.'*P3)\P3.'*d;
%     fd_est(k) = c3(2);
% end
u = ttm(T.core,T.U{1},1);
us = ttm(u,T.U{2},2);
u1 = (double(tenmat(us,3))).';
Psi1 = ((Jn1*u1)\(Jn2*u1)).';
[T1,Phi1] = eig(Psi1);
f2_est = angle(diag(Phi1))/2/pi;
theta_est = asin(-f2_est*2)*180/pi;
Te = rowPermutateMat(M,N);
A = (u1)*conj(T1);
Phi2 = ((Jm1*A)\(Jm2*A)).';
% [T2,Phi2] = eig(Psi2);
f1_est = angle(diag(Phi2))/2/pi;
phi_est = asin(-f1_est*2)*180/pi;
% figure(1)
% hold on
% scatter(phi_est,theta_est)
[phi_est,theta_est] = sort2D(phi_est,theta_est,1);
[phi_est1,theta_est1] = sort2D(phi_est1,theta_est1,1);
[phi1,theta1] = sort2D(phi1,theta1,1);
rmse(nn)  = rmse(nn) + sum(((phi_est.' - phi(1:K)).^2 + (theta_est.' - theta(1:K)).^2))/2/K;
rmse2(nn)  = rmse2(nn) + sum(((phi_est1 - phi(1:K)).^2 + (theta_est1 - theta(1:K)).^2))/2/K;
rmse1(nn) = rmse1(nn) + sum(((phi1.' - phi(1:K)).^2 + (theta1.' - theta(1:K)).^2))/2/K;
end
rmse(nn) = sqrt(rmse(nn)/MM)
rmse1(nn) = sqrt(rmse1(nn)/MM)
rmse2(nn) = sqrt(rmse2(nn)/MM)
end
figure
semilogy(SNR,rmse,'-o',SNR,rmse1,'-s',SNR,rmse2,'-d')
% plot(SNR,mag2db(rmse1),'-o',SNR,mag2db(rmse),'-s',SNR,mag2db(rmse2),'-d')