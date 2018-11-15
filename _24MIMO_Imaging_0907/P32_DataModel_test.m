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
Q = size(targetPos,1);
R_r = c*tc/2;
R_min = R0-R_r/2;
R_max = R0+R_r/2;
T_min = 2*R_min/c;

% targetPos1 = targetPos - ones(I,1)*targetIniPos;
Qi = targetPos.';
targetVel = targetVel.';
P = targetIniPos(:);
qi = Qi - P;
u0 = P/norm(P);
TantPos = [(0:Nt-1)*dt-(Nt-1)*dt/2;zeros(1,Nt);zeros(1,Nt)];  %Transmit ULA and Receive ULA located at x-axis, centered the origin
RantPos = [(0:Nr-1)*dr-(Nr-1)*dr/2;zeros(1,Nr);zeros(1,Nr)];
for nr = 1:Nr
    for nt = 1:Nt
        VantPos(:,(nr-1)*Nt+nt) = (RantPos(:,nr) + TantPos(:,nt))/2;
    end
end
d = VantPos(1,2) - VantPos(1,1);
Nv = Nt*Nr;
for q = 1:Q
    u(q) = qi(:,q).'*u0;
    q0 = qi(:,q) - u(q)*u0;    
    for nv = 1:Nv
        dv = VantPos(:,nv);
        v(q,nv) = q0.'*dv/R0;
        e(q,nv) = dv.'*dv*u(q)/2/R0/R0;
        deltaR(q,nv) =  rangeangle(dv,P)- rangeangle(dv,Qi(:,q));
        deltaR1(q,nv) = v(q,nv) + e(q,nv) - u(q);
%         l(q,nv) = ceil((rangeangle(dv,Qi(:,q))-R_min)/d0);
%         for k = 1:K
%             t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));
%             l1(q,nv,k) = find(rectpuls((t-(k-1)*T-2*rangeangle(dv,Qi(:,q))/c)/T-1/2+1e-8)>0,1);
%         end
%         l1(q,nv,:)
    end
    d(q) = (mean(deltaR(q,2:end) - deltaR(q,1:end-1)))/lambda;
    d1(q) = (mean(deltaR1(q,2:end) - deltaR1(q,1:end-1)))/lambda;
    d2(q) = (mean(v(q,2:end) - v(q,1:end-1)))/lambda;
    d3(q) = (mean(e(q,2:end) - e(q,1:end-1)))/lambda;
end