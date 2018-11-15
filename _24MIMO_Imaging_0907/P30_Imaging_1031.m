clear;clc;close all;

load('data/RxSignal_plane_1106.mat');
SNR = 5;
% Echo = RxSig/sqrt(N);
Echo = awgn(RxSig/sqrt(N),SNR,'measured');
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

NFFT = 1024;
Y1 = reshape(Y(:,1),[N0,Nv]);
for nv = 1:Nv
    X11(:,nv) = fftshift(ifft(Y1(:,nv),NFFT));
end
for nf = 1:NFFT
    X111(nf,:) = fftshift(ifft(X11(nf,:),NFFT));
end

u = (-(NFFT/2-1):NFFT/2)/NFFT*N0*d0;
% theta0 = A0(1)*pi/180;
% phi0 = A0(2)*pi/180;
v = (-(NFFT/2-1):NFFT/2)/NFFT*lambda*R0/2/d;
% x = v - u*u0(1);
% y =  (- x*u0(1)-u)/u0(2);
% x0 = u0*cos(theta0)-v0*sin(theta0);
% y0 = u0*sin(theta0)+v0*cos(theta0);
imagesc(u,v,abs(X111))
axis xy
% xlim([-10,15])
% ylim([-15,10])
% figure
% scatter(u,v)
[phi,psi] = sort2D(phi,psi,1);
u1 = phi*N0*d0;
v1 = psi*lambda/2*R0/d;
figure
scatter(u1,v1,12)
tic
[phi_est,psi_est,s_est] = Unitary_ESPRIT_2D1115(Y,N0,Nv,I);
toc
[phi_est,psi_est] = sort2D(phi_est,psi_est,1);
tic
[phi_est1,psi_est1] = Unitary_ESPRIT_2D(Y,N0,Nv,I);
toc
[phi_est1,psi_est1] = sort2D(phi_est1,psi_est1,1);
u_est = phi_est*N0*d0;
v_est = psi_est*lambda/2*R0/d;
figure
scatter(u_est,v_est,s_est);
x_est = v_est - u_est*u0(1);
y_est =  (- x_est*u0(1)-u_est)/u0(2);
figure
scatter(x_est,y_est,s_est);
u_est1 = phi_est1*N0*d0;
v_est1 = psi_est1*lambda/2*R0/d;
figure
scatter(u_est1,v_est1);
x_est1 = v_est1 - u_est1*u0(1);
y_est1 =  (- x_est1*u0(1)-u_est1)/u0(2);
figure
scatter(x_est1,y_est1);



% tic
%
% toc

    