clear;clc;close all;

c = 3e8;
[Fc,Bw,fs,N,K,df,tb,tc,T,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
lambda = c/Fc;
N_T = fs*T;
N_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
[targetRCS,targetPos,targetVel,targetIniPos] = TargetPara;
[R0,A0] = rangeangle(targetIniPos.',[0;0;0]);
d0 = c/2/Bw;
I = size(targetPos,1);
R_r = c*tc/2;
R_min = R0-R_r/2;
R_max = R0+R_r/2;
T_min = 2*R_min/c;

% targetPos1 = targetPos - ones(I,1)*targetIniPos;
Q = targetPos.';
targetVel = targetVel.';
P = targetIniPos(:);
a = Q - P;
u0 = P/norm(P);
TantPos = [(0:Nt-1)*dt;zeros(1,Nt);zeros(1,Nt)];
RantPos = [zeros(1,Nr);(0:Nr-1)*dr;zeros(1,Nr)];
% TantPos = [(0:Nt-1)*dt-(Nt-1)*dt/2;zeros(1,Nt);zeros(1,Nt)];  %Transmit ULA and Receive ULA located at x-axis, centered the origin
% RantPos = [(0:Nr-1)*dr-(Nr-1)*dr/2;zeros(1,Nr);zeros(1,Nr)];
% for nr = 1:Nr
%     for nt = 1:Nt
%         VantPos(:,(nr-1)*Nt+nt) = (RantPos(:,nr) + TantPos(:,nt))/2;
%     end
% end
% d = VantPos(1,2) - VantPos(1,1);
Nv = Nr*Nt;
P1Q1 = (a - ((a.'*u0)*u0.').').'*[1;0;0];

theta0 = A0(1)*pi/180;
varphi0 = A0(2)*pi/180;
TP = rangeangle(TantPos,P);
PR = rangeangle(RantPos,P);
for q = 1:I
    TQ(q,:) = rangeangle(TantPos,Q(:,q));
    QR(q,:) = rangeangle(RantPos,Q(:,q));
    u(q) = -a(:,q).'*u0;
    v1(q) = (a(:,q)+u(q)*u0).'*[1;0;0];
    v2(q) = (a(:,q)+u(q)*u0).'*[0;1;0];
end
% dR = TP-TQ;
% mean(dR(:,2:end) - dR(:,1:end-1),2)/lambda
psi1 = 2*v1*dt/lambda/R0;
psi2 = 2*v2*dr/lambda/R0;
phi = u/d0/N;

ank = exp(1i*2*pi*ceil(rand(N,K)*4)/4);
RxSig = zeros(K*N_T+N_c,nr);
sigma = rand(1,I)*0.3;
h1 = waitbar(0,'生成数据');
for nr = 1:Nr
    for k = 1:K
        t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));
        for nt = 1:Nt
%             nv = (nr-1)*Nt+nt;
            for q = 1:I
                tau_qtr = (TQ(q,nt)+QR(q,nr))/c;
%                 tau_0(nv) = 2*(R0-(nv-1)*dv*cos(theta_0))/c;
%                 delta_tau(nv) = 2*(u(q)+sin(theta_0)*v(q)/R0*dv*(nv-1))/c;
                sig = targetRCS(q)*exp(1i*2*pi*sigma(q)*(k-1))...
                    *sum(kr(exp(1i*2*pi*Fc*(t-tau_qtr))...
                    .*rectpuls((t-(k-1)*T-tau_qtr)/T-1/2+1e-8),...
                    diag(ank(nt:Nt:end,k))*exp(1i*2*pi*((nt:Nt:N)-1).'...
                    *df*(t-tau_qtr-k*tc))),1);
                RxSig((k-1)*N_T+(1:(N_T+N_c)),nr) = RxSig((k-1)*N_T+(1:(N_T+N_c)),nr)...
                +sig.';
                clear sig;
            end
        end
    end
    waitbar(nr/Nr);
end
delete(h1);
% RxSig = awgn(RxSig,SNR,'measured');
save('data\RxSignal_plane_1207.mat');
% tmp = (dR(:,2:end) - dR(:,1:end-1));%/lambda;
% tmp1 = (dR1(:,2:end) - dR1(:,1:end-1))/lambda;