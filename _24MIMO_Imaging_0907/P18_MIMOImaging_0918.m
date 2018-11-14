clear;clc;close all;

load('RxSignal_03.mat');
N0 = N/Nt;
Nv = Nt*Nr;     %number of virtual element
% X = zeros(Ntb*K,Nr);
% X = [];
RxSig = awgn(RxSig,-10,'measured');
X = zeros(N0*Nv,K);
% F = exp(-1i*2*pi*(0:N-1).'*(0:Ntb-1)/Ntb);

% R_v0 = rangeangle(VantPos,targetIniPos.');
% R_t0 = rangeangle(TantPos,targetIniPos.');
% R_0r = rangeangle(RantPos,targetIniPos.');
tau0 = 2*R0/c;
%% 去CP,解调
for k = 1:K
%     RxSig(:,:,k) = awgn(RxSig(:,:,k),5,'measured');
    t = T_min + (k-1)*T + tc + (0:1/fs:(tb-1/fs));
    t = t(:);
    for nr = 1:Nr
        tmp = exp(-1i*2*pi*Fc*t).*RxSig((k-1)*L_T+((L_c+1):L_T),nr);%RxSig((L_c+1):L_T,nr);   %去CP
        tmpI = fft(tmp)./a(:,k)/N;  %解调，解码
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            tau0(nv) = 2*Rv0(nv)/c;
            tmp0 = tmpI(nt:Nt:end).*exp(1i*2*pi*((nt:Nt:N)-1)*df*(tau0(nv)-T_min)).'...
                *exp(1i*2*pi*Fc*tau0(nv));
            %             plot(abs(ifft(tmp0)))
            X((nv-1)*N0+(1:N0),k)=tmp0;
        end
    end
    
end
X1 = reshape(X(:,1),[N0,Nv]);
NFFT = 1024;
for nv = 1:Nv
    X11(:,nv) = fftshift(ifft(X1(:,nv),NFFT));
end
for n0 = 1:NFFT
    X111(n0,:) = fftshift(fft(X11(n0,:),NFFT));
end
% imagesc(abs(X111))

% tic
%
% toc
uu = c*tb/2/Nt*(-(NFFT/2-1):NFFT/2)/NFFT;
vv = lambda/2/delta_A*(-(NFFT/2-1):NFFT/2)/NFFT;
xx = uu*cosd(A0(1))+vv*sind(A0(1));
yy = uu*sind(A0(1))+fliplr(vv)*cosd(A0(1));
figure
imagesc(xx,yy,abs(X111))
tic
[ku1,kv1,S1] = Unitary_ESPRIT_2D1(X,N0,Nv,Q);
toc
[ku1,inds1] = sort(-ku1);
kv1 = kv1(inds1);
[ku,inds] = sort(ku);
kv = kv(inds);
% S = S/N0;
u_est = ku1/Nt/df*c/2;
v_est = kv1/delta_A*c/Fc/2;
x_est = u_est*cosd(A0(1))+v_est*sind(A0(1));
y_est = u_est*sind(A0(1))-v_est*cosd(A0(1));

ku = ku.';
kv = kv.';
u = ku/Nt/df*c/2;
v = kv/delta_A*c/Fc/2;
x = u*cosd(A0(1))+v*sind(A0(1));
y = u*sind(A0(1))-v*cosd(A0(1));
% tic
% [f11,f21,S1] = ESPRIT_2D(X,N0,Nv,Q);
% toc
% u_est1 = -f11/Nt/df*c/2;
% v_est1 = f21/delta_A*c/Fc/2;
% x_est1 = u_est1*cosd(A0(1))+v_est1*sind(A0(1));
% y_est1 = u_est1*sind(A0(1))-v_est1*cosd(A0(1));
figure
% scatter([0;y_est],[0;x_est],[],[-inf;S],'filled')
% colormap hot
scatter(x_est,y_est,20*log(S1),'filled')
hold on
% % plot(x_est1,y_est1,'s')
plot(x,y,'*')
% x = R_est.*cosd(th_est);
% y = R_est.*sind(th_est);
% [x y zeros(8,1)] - ones(8,1) * targetIniPos;

