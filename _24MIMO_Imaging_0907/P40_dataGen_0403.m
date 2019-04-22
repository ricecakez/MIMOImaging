clear;clc;close all;

c = 3e8;
SNR = 5;
[f0,lambda,B,fs,N,K,df,tb,tc,ts,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
deltaR = c/2/B;
L_T = fs*ts;
L_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
Nv = Nt*Nr;
[rho,SP,V,P] = TargetPara;
I = size(SP,1);
Q = SP.';
P = P.';
q = Q - P;
T = [dt*(0:Nt-1)-dt*(Nt-1)/2;zeros(2,Nt)];
R = [dr*(0:Nr-1)-dr*(Nr-1)/2;zeros(2,Nr)];
OP = rangeangle(P);
Rc = c*tc/2;
Rmin = OP-Rc/2;
Tmin = 2*Rmin/c;
TP = rangeangle(T,P);
PR = rangeangle(R,P);
OQ = rangeangle(Q);
u0 = P/OP;
ui = q.'*u0;
vi = (q-u0*ui.').'*[1;0;0];
for nt = 1:Nt
    for nr = 1:Nr
        nv = (nr-1)*Nt+nt;
%         V = [(nv-1)*dt/2;0;0];
%         R_VP(nv,1) = 2*rangeangle(V,P);
%         R_VQ(nv,:) = 2*rangeangle(SP,V);
%         dR(nv,:) = R_VP(nv)-R_VQ(nv,:) - ui.'*(((nt-1)*dt)^2+((nr-1)*dr)^2)/2/R0^2;
        TPR(nv,1) = TP(nt) + PR(nr);
        TQR(nv,:) = rangeangle(Q,T(:,nt))+rangeangle(Q,R(:,nr));
        dR0(nv,:) = TPR(nv) - TQR(nv,:);
        dR1(nv,:) = -2*ui+vi*(T(1,nt)+R(1,nr))/OP...
            +ui*(T(1,nt)^2+R(1,nr)^2)/2/OP^2;
    end
end
for i = 1:I
    dd0(:,i) = dR0(2:end,i) - dR0(1:end-1,i);
    dd1(:,i) = dR1(2:end,i) - dR1(1:end-1,i);
    figure
    plot(1:Nv-1,dd0(:,i),1:Nv-1,dd1(:,i));
end
tau_nv = (TPR-dR1)/c;
tau_nv1 = TQR/c;
a = exp(1i*2*pi*ceil(rand(N,K)*4)/4);   %QPSK
B = diag(rho)*exp(1i*randn(I,K));
Y = zeros(N*Nr,K);
eps = 1e-8;
h1 = waitbar(0,'生成数据');
for k = 1:K
    t = Tmin + (k-1)*ts + tc + (0:1/fs:(tb-1/fs));
    for nr = 1:Nr
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
%             tau_nv = R_TQR(nv,:)/c;
            m = nt:Nt:N;
            for i = 1:I
                Y((nr-1)*N+(1:N),k) = Y((nr-1)*N+(1:N),k)...
                    + B(i,k)*exp(1i*2*pi*f0*(t.'-tau_nv(nv,i)))...
                    .*sum(diag(a(m,k))*exp(1i*2*pi*(m-1).'*df*(t-tau_nv(nv,i)-k*tc)),1).'...
                    .*rectpuls((t.'-(k-1)*ts-tau_nv(i))/ts-1/2+eps);
            end
%             nv = 
        end
    end
     waitbar(k/K);
end
delete(h1);
% RxSig = awgn(RxSig,SNR,'measured');
save('data\RxSignal_targetPlane_1_0409.mat');

