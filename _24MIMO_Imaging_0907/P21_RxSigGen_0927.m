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
R_min = R0-R_r/2;
R_max = R0+R_r/2;
T_min = 2*R_min/c;
a = exp(1i*2*pi*ceil(rand(N,K)*4)/4);   %QPSK
L_T = round(T*fs);
Sig_nv = zeros(L_T*K+L_c,1);
RxSig = zeros(L_T,Nr,K);
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
    u(q) = targetPos1(q,1)*cosd(A0(1))+targetPos1(q,2)*sind(A0(1));
    v(q) = targetPos1(q,1)*sind(A0(1))-targetPos1(q,2)*cosd(A0(1));
    ku(q) = 2*Nt*u(q)/c/tb;
    kv(q) = 2*v(q)*delta_A/lambda;
end
h1 = waitbar(0,'生成数据');
Rv0 = rangeangle(VantPos,targetIniPos.');
sigma = rand(1,Q)*0.3;
for nr = 1:Nr
    for k = 1:K
        t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            for q = 1:Q
                Rv(nv,q) = Rv0(nv) + u(q) - (nv-1)*delta_A*v(q);
                tau_vq(nv,q) = 2*Rv(nv,q)/c;
                sig = targetRCS(q)*exp(1i*2*pi*sigma(q)*(k-1))...
                    *sum(kr(exp(1i*2*pi*Fc*(t-tau_vq(nv,q)))...
                    .*rectpuls((t-(k-1)*T-tau_vq(nv,q))/T-1/2+1e-8),...
                    diag(a(nt:Nt:end,k))*exp(1i*2*pi*((nt:Nt:N)-1).'...
                    *df*(t-tau_vq(nv,q)-k*tc))));
                Sig_nv((k-1)*L_T+(1:(L_T+L_c))) = Sig_nv((k-1)*L_T+(1:(L_T+L_c)))...
                    +sig.';
                clear sig;
            end
        end
        RxSig(:,nr,k) = Sig_nv((k-1)*L_T+(1:L_T));
    end
    waitbar(nr/Nr);
end
delete(h1);
% RxSig = awgn(RxSig,SNR,'measured');
save('RxSignal_plane4.mat');
