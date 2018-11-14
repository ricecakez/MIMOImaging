clear;clc;close all;

% [ L0,Omega0,dOmega,V0,a0] = ParametersTarget();
% [ F0,F_sample,B,N,df,tb,tc,T,PRF,K,Tp,c] = ParametersSystem();
load ('ReturnSignal.mat');
figure
imagesc(abs(signal_return)/max(max(abs(signal_return)))),colormap(gray);
title('回波信号图');
dR = c/2/B;
dV = c/F0/2/T/K;
L_range = c*tc/2;
L_min = L0-L_range/2;     %最早接受数据对应的目标距离，也即最小的接收数据的距离

%% 去CP

ic = round(tc*F_sample);
iT = round(T*F_sample);
iR = size(signal_return,1);
signal_process = signal_return(ic+1:iT,:);
t_rec = 2*L_min/c+(ic+1:iT).'/F_sample;
for k = 1:K
    tmp = fft(signal_process(:,k).*exp(-1i*2*pi*F0*t_rec));
    y(:,k) = tmp(1:N);
    y1(:,k) = y(:,k)./a(:,k);
    Y1(:,k) = ifft(y1(:,k));
end
figure,imagesc(0:K-1,L_min+dR*(0:N-1),abs(Y1)/max(max(abs(Y1)))),colormap(gray);
title('距离压缩结果');
y11 = reshape(y1,[N*K,1]);
[f1,f2] = Unitary_ESPRIT_2D1(y11,N,K,8);
tau_est = -f1/2/df;
R_est = tau_est*c/2 - L0 + L_min;
f2*c/4/Omega*T
% R_est - L0 + L_min
% t = linspace(-0.5,0.5,512);
% f = linspace(-1e-2,1e-2,512);
% [t_e,f_e] = RD_MUSIC(r,p,t,f)
% s = MUSIC_est2D(y1,[8,40],t,f,40,40);
for n = 1:N
    Y(n,:) = ifft(Y1(n,:),1024);
end
output = abs(Y)/max(max(abs(Y)));
figure,contour(output),colormap(jet);  %contour在图像较为复杂的时候计算时间非常长，当图像复杂时慎用
title('成像结果轮廓图');
figure,imagesc((0:1024-1)*dV*K/1024/Omega0,L_min+dR*(0:N-1),output),colormap(gray),colormap(hot);
title('成像结果黑白图');
xlim([-10,10])

L = 20;
% f = linspace(-1e-2,1e-2,512);
% figure
for n = 250:261
    [f(n,:) G(n,:)]= Unitary_ESPRIT(y1(n,:),40,15);
    [f(n,:) inds] = sort(f(n,:));
    G(n,:) = G(n,inds);
%     output1(n,:) = S;
%     hold on
%     plot(f*c/F0/T/Omega0/2,output1(n,:))
%     Ryy = Y1(n,:)'*Y1(n,:);
%     RSM = spsmooth(Ryy,L,'fb');
%     doasm(n,:) = rootmusicdoa(RSM,18);
end
[i,j] = find(G==max(max(G,15)));
for p = 1:15
    ts(p) = i(p)-1;
    fs(p) = f(i(p),j(p));
    Gs(p) = G(i(p),j(p));
end
% f(i,j)
figure
plot(L_min+dR*ts,fs*c/F0/T/Omega0/2,'o')
imagesc(L_min+dR*(0:N-1)-L0,f*c/F0/T/Omega0/2,output1.'),colormap(hot);
L_max = tc*c/4;
ylim([-10 10])
xlim([-L_max L_max])



% p1 = 50;
% p2 = 50;
% l = 0;
% Rxx = zeros(p1*p2);
% J = fliplr(eye(p1*p2));
% for n = 1:N-p1+1
%     for k = 1:K-p2+1   
%         l = l+1;
%         tmp = reshape(y1(n:n+p1-1,k:k+p2-1),[p1*p2,1]);
%         Rll = tmp*tmp';
%         Rxx = Rll + J*conj(Rxx)*J;
%     end
% end
% Rxx = Rxx/2/l;
% y2 = reshape(y1,[N*K,1]);
% R = y2*y2';
% L = 100;
% RSM = spsmooth(R,L);
% doasm = rootmusicdoa(RSM,18)