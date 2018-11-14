clear;clc;close all;
c = 3e8;
SNR = 5;
[Fc,Bw,fs,N,K,df,tb,tc,T,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
L_T = fs*T;
L_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
[targetRCS,targetPos,targetVel,targetIniPos] = TargetPara;
Q = size(targetPos,1);

targetPos1 = targetPos - ones(Q,1)*targetIniPos;
targetPos = targetPos.';
targetVel = targetVel.';

Q = 2;
TantPos = [(0:Nt-1)*dt;zeros(1,Nt);zeros(1,Nt)];  %Transmit ULA and Receive ULA located at x-axis, centered the origin
RantPos = [(0:Nr-1)*dr;zeros(1,Nr);zeros(1,Nr)];
R_r = c*tc/2;
[R0,A0] = rangeangle(targetIniPos.',zeros(3,1));
R_min = R0-R_r/2;
R_max = R0+R_r/2;
T_min = 2*R_min/c;
a = exp(1i*2*pi*ceil(rand(N,K)*4)/4);   %QPSK
L_T = round(T*fs);
RxSig = zeros(L_T*K+L_c,Nr);

t = T_min + (0:1/fs:(T+tc-1/fs));
for nr = 1:Nr
    for nt = 1:Nt
        VantPos(:,(nr-1)*Nt+nt) = (RantPos(:,nr) + TantPos(:,nt))/2;
    end
end
Nv = Nr*Nt;
dv = dt/2;
Lv = Nv*dv;
delta_A = acos((2*R0^2-Lv^2)/2/R0^2)/Nv;
% [R_v0,A_v0] = rangeangle(VantPos,targetIniPos.');
% A_v0(1,:) = 180+A_v0(1,:);
% delta_A = mean(A_v0(1,2:end) - A_v0(1,1:end-1))*pi/180;
% delta_A0 = (A_v0(1,end) - A_v0(1,1))/Nr/Nt*pi/180;
% for q = 1:Q
%     [R_rq, A_rq] = rangeangle(RantPos,targetPos(:,q));
%     [R_tq, A_tq] = rangeangle(TantPos,targetPos(:,q));
%     [R_vq, A_vq] = rangeangle(VantPos,targetPos(:,q));
%     R_vq = R_vq * 2;
%     for nr = 1:Nr
%         for nt = 1:Nt
%             R_tqr((nr-1)*Nt+nt) = R_rq(nr)+R_tq(nt);
%         end
%     end
%     d1(:,q) = R_tqr(2:end) - R_tqr(1:end-1);
%     d2(:,q) = R_vq(2:end) - R_vq(1:end-1);
%     f2(q) = 2*Fc*mean(d1(:,q))/c;
% end
h1 = waitbar(0,'生成数据');
for k = 1:K
    for q = 1:Q
                u(q) = targetPos1(q,1)*cosd(A0(1))+targetPos1(q,2)*sind(A0(1));
                v(q) = targetPos1(q,1)*sind(A0(1))-targetPos1(q,2)*cosd(A0(1));
                ku(q) = 2*Nt*u(q)/c/tb;
                kv(q) = 2*v(q)*delta_A*Fc/c;
        %         r = R_v0 + u*cos((0:Nr*Nt-1)*delta_A) - v*sin((0:Nr*Nt-1)*delta_A);
        %         dtr = -dt*sind(A0(1));
        %         dR_v0 = R_v0(2:end)-R_v0(1:end-1);
        %         Rv = 2*rangeangle(VantPos,targetPos(:,q));
        %         dr = 2*(mean(r(2:end) - r(1:end-1)));
        [R_rq, A_rq] = rangeangle(RantPos,targetPos(:,q));
        [R_tq, A_tq] = rangeangle(TantPos,targetPos(:,q));
        %         for nr = 1:Nr
        %             for nt = 1:Nt
        %                 R_tqr((nr-1)*Nt+nt) = R_rq(nr)+R_tq(nt);
        %             end
        %         end
        %         dRtqr = R_tqr(2:end) - R_tqr(1:end-1);
        for nr = 1:Nr
            for nt = 1:Nt
                tau_tqr = (R_tq(nt) + R_rq(nr))/c;
%                 tau_tqr1 = 2*R_tqr((nr-1)*Nr+nt)/c;
%                 Lntnr((nr-1)*Nt+nt,q) = -2*Fc*(tau_tqr-2*R0/c);
%                 Ln(nt,nr,q) = ceil((T_min-tau_tqr)*fs);
                sig = targetRCS(q)*sum(kr(diag(a(nt:Nt:end,k))*exp(1i*2*pi*((nt:Nt:N)-1).'*df*(t-tau_tqr-tc)),...
                    exp(1i*2*pi*Fc*(t-tau_tqr)).*rectpuls((t-tau_tqr)/T-1/2+1e-8)));
                RxSig((k-1)*L_T+(1:(L_T+L_c)),nr) = RxSig((k-1)*L_T+(1:(L_T+L_c)),nr)...
                    +sig.';
                clear sig;
            end
        end
    end
    targetPos = targetPos + targetVel*T;
    waitbar(k/K);
end
delete(h1);
RxSig = awgn(RxSig,SNR,'measured');
save('RxSignal.mat');
