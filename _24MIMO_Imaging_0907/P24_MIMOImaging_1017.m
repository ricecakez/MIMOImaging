clear;clc;close all;

load('RxSignal_plane_new.mat');
% RxSig = awgn(RxSig,-10,'measured');

%% Preprocessing
for k = 1:K
    t = T_min + (k-1)*T + tc + (0:1/fs:tb-1/fs);
    t = t(:);
%     Xk = awgn(RxSig(:,:,k),SNR,'measured');
    for nr = 1:Nr
        y = RxSig((k-1)*L_T+((L_c+1):L_T),nr).*exp(-1i*2*pi*Fc*t);
        y1 = fft(y)./a(:,k)/N;
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            tau0 = 2*(R0-(nv-1)*dv*cos(theta_0))/c;
            y2 = y1(nt:Nt:end).*exp(1i*2*pi*((nt:Nt:N)-1)*df*(tau0-T_min)).'...
                *exp(1i*2*pi*Fc*tau0);            
            Y((nv-1)*N0+(1:N0),k) = y2;
        end
    end
end
Y1 = reshape(Y(:,1),[N0,Nv]);
NFFT = 1024;
for nv = 1:Nv
    X11(:,nv) = fftshift(ifft(Y1(:,nv),NFFT));
end
for n0 = 1:NFFT
    X111(n0,:) = fftshift(ifft(X11(n0,:),NFFT));
end
uu = c*tb/2/Nt*(-(NFFT/2-1):NFFT/2)/NFFT;
vv = lambda/2/delta_A*(-(NFFT/2-1):NFFT/2)/NFFT;
xx = uu*cosd(A0(1))+vv*sind(A0(1));
yy = uu*sind(A0(1))-fliplr(vv)*cosd(A0(1));
figure
imagesc(uu,vv,abs(X111))
[ku1,kv1,S1] = Unitary_ESPRIT_2D1(Y,N0,Nv,Q);
[ku1,inds1] = sort(ku1);
[ku,inds] = sort(phi.');
kv = psi(inds).';
u = u(inds).';
v = v(inds).';
kv1 = kv1(inds1);
u1 = ku1/Nt/df*c/2;
v1 = kv1/dv/sin(theta_0)*lambda*R0/2;
x1 = u1*cos(theta_0)+v1*sin(theta_0);
y1 = -u1*sin(theta_0)+v1*cos(theta_0);
figure
scatter(x1,y1)
