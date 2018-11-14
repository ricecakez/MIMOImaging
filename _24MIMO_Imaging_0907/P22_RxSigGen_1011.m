clear;clc;close all;
c = 3e8;
SNR = -5;
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
a = exp(1i*2*pi*ceil(rand(1,K)*4)/4);   %QPSK
L_T = round(T*fs);
RxSig = zeros(L_T*K+L_c,Nr);
for nt = 1:Nt
    for nr = 1:Nr
        VantPos(:,(nr-1)*Nt+nt) = (RantPos(:,nr) + TantPos(:,nt))/2;
    end
end
Nv = Nr*Nt;
dv = VantPos(1,2)-VantPos(1,1);
Lv = Nv*dv;
delta_A = dv/R0;%acos((2*R0^2-Lv^2)/2/R0^2)/Nv;
for q = 1:Q
    u(q) = targetPos1(q,1)*cosd(A0(1))+targetPos1(q,2)*sind(A0(1));
    v(q) = targetPos1(q,1)*sind(A0(1))-targetPos1(q,2)*cosd(A0(1));
    %         kv1(q) = 2*dv*cosd(A0(1))*v(q)/(u(q)+R0)/lambda;
    ku(q) = 2*Nt*u(q)/c/tb;
    kv(q) = 2*v(q)*delta_A/lambda;
end

h1 = waitbar(0,'生成数据');
Rv0 = rangeangle(VantPos,targetIniPos.');
sigma = rand(1,Q)*0.3;

for k = 1:K
    t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));
    for q = 1:Q
        %         Rt(:,q) = rangeangle(TantPos,targetPos(:,q));
        %         Rr(:,q) = rangeangle(RantPos,targetPos(:,q));
        
        for nv = 1:Nv
            nr = floor((nv-1e-8)/Nt)+1;
            nt = mod(nv,Nt);
            if nt == 0
                nt = Nt;
            end
            Rv(nv,q) = Rv0(nv)+u(q)-(nv-1)*delta_A*v(q);
            tau_vq(nv,q) = 2*Rv(nv,q)/c;
            %             for nr = 1:Nr
            %                 R_q((nt-1)*Nr+nr,q) = Rt(nt,q) + Rr(nr,q);
            %                 tau_q((nt-1)*Nr+nr,q) = R_q((nt-1)*Nr+nr,q)/c;
            
            %                 Rv((nt-1)*Nr+nr,q) = Rv0((nt-1)*Nr+nr) + u(q) - ((nt-1)*Nr+nr-1)*delta_A*v(q);
            %                 tau_vq((nt-1)*Nr+nr,q) = 2*Rv((nt-1)*Nr+nr,q)/c;
            %                 tau_vq0((nt-1)*Nr+nr,q) = tau_vq((nt-1)*Nr+nr,q) - 2*Rv0((nt-1)*Nr+nr)/c;
            sig = targetRCS(q)*exp(1i*2*pi*sigma(q)*(k-1))...
                *sum(kr(exp(1i*2*pi*Fc*(t-tau_vq(nv,q)))...
                .*rectpuls((t-(k-1)*T-tau_vq(nv,q))/T-1/2+1e-8),...
                diag(a(:,k))*exp(1i*2*pi*((nt:Nt:N)-1).'...
                *df*(t-tau_vq(nv,q)-k*tc))));
            RxSig((k-1)*L_T+(1:(L_T+L_c)),nr) = RxSig((k-1)*L_T+(1:(L_T+L_c)),nr)...
                +sig.';
            clear sig;
        end
        
    end
    
    %     end
    
    %     targetPos = targetPos + targetVel*T;
    waitbar(k/K);
end
delete(h1);
save('RxSignal_04.mat');
% RxSig = awgn(RxSig,SNR,'measured');

