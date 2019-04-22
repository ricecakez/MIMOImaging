clear;clc;close all;

load('data\RxSignal_plane_1226.mat');
SNR = 15;
ms = sort_matrix([phi,psi,phi+psi],'ascend',1);
phi = ms(:,1);
psi = ms(:,2);
u = phi*N*d0;
v = psi*lambda*R0/d/2;
x = v - u*u0(1);
y =  (- x*u0(1)-u)/u0(2);
Echo = RxSig/sqrt(N);
% Echo = awgn(RxSig/sqrt(N),SNR,'measured');
for k = 1:K
    t = T_min + (k-1)*T + (tc:1/fs:(T-1/fs));
    t = t(:);
    for nr = 1:Nr
        yc = exp(-1i*2*pi*Fc*t).*Echo((k-1)*N_T+((N_c+1):N_T),nr);
        tmp = fft(yc)./ank(:,k)/sqrt(N);
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            TP = rangeangle(TantPos(:,nt),P);
            PR = rangeangle(P,RantPos(:,nr));
            Phi0 = exp(1i*2*pi*((nt:Nt:N)-1).'*df*((TP+PR)/c-T_min));
            Psi0 = exp(1i*2*pi*(TP+PR)/lambda);
            Y((nv-1)*N0+(1:N0),k) = tmp(nt:Nt:N).*Phi0*Psi0;
        end
    end
end

NFFT = 64;
Y1 = reshape(Y(:,1),[N0,Nv]);
for nv = 1:Nv
    X11(:,nv) = fftshift(fft(Y1(:,nv),NFFT))/Nv;
end
for nf = 1:NFFT
    X111(nf,:) = fftshift(fft(X11(nf,:),NFFT))/N0;
end
u_ = (-NFFT/2:(NFFT/2-1))*d0*N0/NFFT;
v_ = (-NFFT/2:(NFFT/2-1))/Nv*lambda*R0/2/d*Nv/NFFT;

figure
imagesc(v_,u_,abs(X111))
axis xy
xlabel('$v/{\mathrm{m}}$','Interpreter','latex')
ylabel('$u/{\mathrm{m}}$','Interpreter','latex')
colormap('hot')

[m,n,S_1] = find_peak_2D(abs(X111),I);
u_1 = u_(m);
v_1 = v_(n);
m_s = sort_matrix([u_1.',u_1.'+v_1.',v_1.',S_1.'],'ascend',1);
u_1 = m_s(:,1);
v_1 = m_s(:,3);
S_1 = m_s(:,4);
x_1 = v_1 - u_1*u0(1);
y_1 = (- x_1*u0(1)-u_1)/u0(2);
figure
scatter(v_1,u_1,mag2db(S_1),'k','filled')
figure
scatter(x_1,y_1,mag2db(S_1),'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
[f1,f2,S1] = Uni_Succ_ESPRIT(Y,N0,Nt,Nr,I);
m_s = sort_matrix([f1,f1+f2,f2,S1],'ascend',1);
f1 = m_s(:,1);
f2 = m_s(:,3);
S1 = m_s(:,4);
u1 = f1*N*d0;
v1 = f2*lambda*R0/d/2;
x1 = v1 - u1*u0(1);
y1 =  (- x1*u0(1)-u1)/u0(2);
figure
% hold on
scatter(x1,y1,S1,'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')

figure
scatter(x,y,mag2db(targetRCS),'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
figure
scatter(v,u,mag2db(targetRCS),'k','filled')
grid on
xlabel('$v/{\mathrm{m}}$','Interpreter','latex')
ylabel('$u/{\mathrm{m}}$','Interpreter','latex')

