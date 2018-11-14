clear;clc;close all;

theta = [10 20 30];
M = 8;
K = 3;
L = 200;
A = exp(-1i*pi*(0:M-1).'*sin(theta*pi/180));
data = randn(K,L);
data = sign(data);
pw = [1 0.8 0.7]';
S = diag(pw)*data;
X = A*S;
nv = ones(1,M);
snr = 10;
X = X + diag(sqrt(nv/snr/2))*(randn(M,L)+1i*randn(M,L));
R = X*X'/L;
% R1 = zeros(M);
% for l = 1:L
%     R1 = R1 + X(:,l)*X(:,l)';
% end
% R1 = R1/L;
X1 = X(1:K,:);
X2 = X(K+1:end,:);
G = R(:,1:K);
H = R(:,K+1:end);
P0 = X2*X1'/(X1*X1');
P = (G'*G)\G'*H;
Q = [P',-eye(M-K)];
Q0 = Q*(Q'*Q)^(-1/2);
[V1,D1] = eig(R);
[d1,ind1] = sort(diag(D1),'descend');
UN1 = V1(:,ind1(K+1:end));
dt = 0.1;
the1 = -90:dt:(90-dt);

Pc = ((G'*G)\G'*H)';
P1 = [eye(K); Pc];
Pa = P1(1:end-1,:);
Pb = P1(2:end,:);
Psir = (Pa'*Pa)\Pa'*Pb;
[V,D] = eig(Psir);
the_est = -asin(angle(diag(D))/pi)*180/pi
for i = 1:length(the1)
    a = exp(-1i*pi*(0:M-1).'*sin(the1(i)*pi/180));
    SP(i) = 1/(a'*(Q'*Q)*a);
    SP1(i) = 1/(a'*(Q0'*Q0)*a);
    MP(i) = (a'*a)/(a'*(UN1*UN1')*a);
end
plot(the1,mag2db(abs(SP)/max(abs(SP))),the1,mag2db(abs(SP1)/max(abs(SP1))),the1,mag2db(abs(MP)/max(abs(MP))))
xlim([0 60])

the_est1 = the1(find_peak(SP,K))
the_est2 = the1(find_peak(SP1,K))
the_est3 = the1(find_peak(MP,K))