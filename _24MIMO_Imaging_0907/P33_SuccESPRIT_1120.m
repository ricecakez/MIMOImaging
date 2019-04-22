clear;clc;close all;

load('data/RxSignal_plane_1106.mat');
f = phi/Nt+psi;
ms = sort_matrix([phi/Nt,psi,f],'ascend',1);
f = ms(:,3);
phi = ms(:,1);
psi = ms(:,2);
u = phi*N*d0;
v = psi*lambda/2*R0/d;
figure
scatter(u,v,12,'ro')
hold on
SNR = 5;
Echo = RxSig/sqrt(N);

% Echo = awgn(RxSig/sqrt(N),SNR,'measured');
for k = 1:K
    t = T_min + (k-1)*T + (tc:1/fs:(T-1/fs));
    t = t(:);
    for nr = 1:Nr
        y = exp(-1i*2*pi*Fc*t).*Echo((k-1)*N_T+((N_c+1):N_T),nr);
        tmp = fft(y)./ank(:,k)/sqrt(N);
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            Phi0 = exp(1i*2*pi*((nt:Nt:N)-1).'*df*(2*VP(nv)/c-T_min));
            Psi0 = exp(1i*2*pi*2*VP(nv)/lambda);
            Y((nv-1)*N0+(1:N0),k) = tmp(nt:Nt:N).*Phi0*Psi0;
        end
    end
end

J3 = kron(kron(eye(Nr),eye(Nt)),selectMat(N0,1));
J4 = kron(kron(eye(Nr),eye(Nt)),selectMat(N0,2));
J5 = kron(selectMat(Nr*Nt,1),eye(N0));
J6 = kron(selectMat(Nr*Nt,2),eye(N0));
% % J5 = kron(eye(Nr),selectMat(N0*Nt,1));
% % J6 = kron(eye(Nr),selectMat(N0*Nt,2));
J1 = kron(kron(eye(Nr),selectMat(Nt,1)),eye(N0));
J2 = kron(kron(eye(Nr),selectMat(Nt,2)),eye(N0));
% J5 = kron(kron(selectMat(Nr,1),eye(Nt)),eye(N0));
% J6 = kron(kron(selectMat(Nr,2),eye(Nt)),eye(N0));
Z = Unitary_transform(Y);
[U,D,~] = svd(Z);
% max(max(U*D*U'-Y*Y'/K))
% d = diag(D);
Es = UniMat(size(Z,1))*U(:,1:I);%*D(1:I,1:I);
Psi_mu = pinv(J3*Es)*(J4*Es);
[T,Phi_mu] = eig(Psi_mu);
rank(T)
[mu,inds] = sort(angle(diag(Phi_mu))/2/pi/Nt);
T = T(:,inds);
A = Es*T;

% C1 = kron(kron(ones(Nr,1),exp(-1i*2*pi*(0:Nt-1).'*mu.')),ones(N0,1));
% A1 = C1.*A;
Phi_mu = pinv(J3*A)*(J4*A);
mu = angle(diag(Phi_mu))/2/pi/Nt;
% [ms,inds] = sort_matrix([mu,nu],'ascend',1);
% mu = ms(:,1);
% nu = ms(:,2);
C = kron(kron(ones(Nr,1),exp(-1i*2*pi*(0:Nt-1).'*mu.')),ones(N0,1));
A1 = C.*A;

nu = f1 - mu;
Phi_nu = pinv(J5*A1)*(J6*A1);
nu1 = angle(diag(Phi_nu))/2/pi;
u_est = mu*N*d0;
v_est = nu1*lambda/2*R0/d;
figure(1)
scatter(u_est,v_est,12,'kx');
% 
% Phi_nu = pinv(J5*A)*(J6*A);
% nu = angle(diag(Phi_nu))/2/pi/Nt;
% nu1 = f12-mu;
% u_est = mu*N*d0;
% v_est = nu1*lambda/2*R0/d;

% figure(1)
% scatter(u_est,v_est,12,'kx');