clear;clc;close all;
c = 3e8;
SNR = 5;
[Fc,Bw,fs,N,K,df,tb,tc,T,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
lambda = c/Fc;
L_T = fs*T;
L_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
[targetRCS,targetPos,targetVel,targetIniPos] = TargetPara;
Q = size(targetPos,1);

targetPos1 = targetPos - ones(Q,1)*targetIniPos;
targetPos = targetPos.';
targetVel = targetVel.';

% Q = 3;
TantPos = [(0:Nt-1)*dt;zeros(1,Nt);zeros(1,Nt)];  %Transmit ULA and Receive ULA located at x-axis, centered the origin
RantPos = [(0:Nr-1)*dr;zeros(1,Nr);zeros(1,Nr)];
R_r = c*tc/2;
[R0,A0] = rangeangle(targetIniPos.',zeros(3,1));
theta_0 = A0(1)*pi/180;
R_min = R0-R_r/2;
R_max = R0+R_r/2;
T_min = 2*R_min/c;
a = exp(1i*2*pi*ceil(rand(N,K)*4)/4);   %QPSK
L_T = round(T*fs);
% Sig_nv = zeros(L_T*K+L_c,1);
RxSig = zeros(L_T*K+L_c,Nr);
for nr = 1:Nr
    for nt = 1:Nt
        VantPos(:,(nr-1)*Nt+nt) = (RantPos(:,nr) + TantPos(:,nt))/2;
    end
end
Nv = Nr*Nt;
dv = dt/2;
Lv = Nv*dv;
delta_A = dv/R0;
for q = 1:Q
    u(q) = targetPos1(q,1)*cos(theta_0)+targetPos1(q,2)*sin(theta_0);
    v(q) = -targetPos1(q,1)*sin(theta_0)+targetPos1(q,2)*cos(theta_0);
    theta_i(q) = theta_0 + v(q)/R0;
    tmp(q) = cos(theta_i(q)) - (cos(theta_0) - v(q)/R0*sin(theta_0))
    phi(q) = -2*Nt/c/tb*u(q);
    psi(q) = 2*dv*sin(theta_0)/lambda/R0*v(q);
%     kv(q) = - cos(A0)*v(q)/(R0+u(q));
end

figure
scatter(targetPos1(:,1),targetPos1(:,2),targetRCS)
figure
scatter(u,v,targetRCS.')
h1 = waitbar(0,'生成数据');
Rv0 = rangeangle(VantPos,targetIniPos.');
sigma = rand(1,Q)*0.3;
for nr = 1:Nr
    for k = 1:K
        t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            Rv(nv,:) = rangeangle(targetPos,VantPos(:,nv));
            for q = 1:Q
                Rv1(nv,q) = R0+u(q)-(cos(theta_0)-sin(theta_0)*v(q)/R0)*dv*(nv-1);
%                 Rv(nv,q) = Rv0(nv) + u(q) - (nv-1)*delta_A*v(q);
                tau_vq(nv,q) = 2*Rv1(nv,q)/c;
                tau_0(nv) = 2*(R0-(nv-1)*dv*cos(theta_0))/c;
                delta_tau(nv) = 2*(u(q)+sin(theta_0)*v(q)/R0*dv*(nv-1))/c;
                sig = targetRCS(q)*exp(1i*2*pi*sigma(q)*(k-1))...
                    *sum(kr(exp(1i*2*pi*Fc*(t-tau_vq(nv,q)))...
                    .*rectpuls((t-(k-1)*T-tau_vq(nv,q))/T-1/2+1e-8),...
                    diag(a(nt:Nt:end,k))*exp(1i*2*pi*((nt:Nt:N)-1).'...
                    *df*(t-tau_vq(nv,q)-k*tc))),1);
                RxSig((k-1)*L_T+(1:(L_T+L_c)),nr) = RxSig((k-1)*L_T+(1:(L_T+L_c)),nr)...
                +sig.';
                clear sig;
            end
        end
    end
    waitbar(nr/Nr);
end
delete(h1);
% RxSig = awgn(RxSig,SNR,'measured');
save('data\RxSignal_plane_1114.mat');
