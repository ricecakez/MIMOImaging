clear;clc;close all;

SNR = 5;
[ f0,df,N,B,fs,tb,tc,T,K,Tp,c,PRI] = OFDMRadarPara;
M = 100;
lambda = c/f0;
L_T = fs*T;
L_c = fs*tc;
[rcs,Tq,vel,To] = TargetParaISAR;
Q = size(Tq,2);
[R0,A0] = rangeangle(To);
[v0,A1] = rangeangle(vel);
A0 = A0(1,1)*pi/180;
A1 = A1(1,1)*pi/180;
A = A0 - A1;
vr = v0*cos(A);
omega = v0*sin(A)/R0;
sig_k = zeros((K-1)*L_T+L_c,1);
X = zeros(L_T,K);
a = exp(1i*2*pi*ceil(rand(N,1)*4)/4);   %QPSK
R_r = c*tc/2;
R_min = R0 - R_r/2;
R_max = R0 + R_r/2;
T_min = 2*R_min/c;

for q = 1:Q
    u(q) = Tq(1,q)*cos(A)+Tq(2,q)*sin(A);
    v(q) = Tq(1,q)*sin(A)-Tq(2,q)*cos(A);
end
scatter(u,v)
for k = 1:K
    t = T_min + (k-1)*T + (0:1/fs:(T+tc-1/fs));    
    Rt = R0 + vr*t;
    for q = 1:Q
        Rq = Rt + u(q) + v(q)*omega*t;
        tau = 2*Rq/c;
        sig = rcs(q)*sum(kr(exp(1i*2*pi*f0*(t-tau))...
            .*rectpuls((t-(k-1)*T-tau)/T-1/2+1e-8),...
            diag(a)*exp(1i*2*pi*(0:N-1).'...
            *df*(t-tau-k*tc))));
        sig_k((k-1)*L_T+(1:(L_T+L_c))) = sig_k((k-1)*L_T+(1:(L_T+L_c)))...
            +sig.';
        clear sig;
    end
    X(:,k) = sig_k((k-1)*L_T+(1:L_T));
end
save('OFDMISAREcho_mat');
