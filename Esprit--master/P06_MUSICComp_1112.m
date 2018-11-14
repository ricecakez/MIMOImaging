clear;clc;close all;

theta = [10,20,30];
M = 8;
L = 200;
A = exp(-1i*pi*(0:M-1).'*sin(theta/180*pi));
K = 3;
% S = randn(K,L);
S0 = exp(1i*2*pi*rand(K-1,L));
S = [S0;S0(1,:)];
P = S*S'/L;
X = A*S;
SNR = 10;
X = awgn(X,SNR,'measured');

dt = 0.1;
theta1 = -90:dt:(90-dt);

R = X*X'/L;
% R1 = A*P*A';
[V1,D1] = eig(R);
[d1,ind1] = sort(diag(D1),'descend');
UN1 = V1(:,ind1(K+1:end));

% rank(R)
Iv = IEM(M);
Rx = R + Iv*conj(R)*Iv;
[V2,D2] = eig(Rx);
[d2,ind2] = sort(diag(D2),'descend');
UN2 = V2(:,ind2(K+1:end));
RA2 = V2(:,ind2(1:K))/diag(d2(ind2(1:K)))*V2(:,ind2(1:K))';

Kelm = 6;
RS = ssp(R,Kelm);
[V3,D3] = eig(RS);
[d3,ind3] = sort(diag(D3),'descend');
UN3 = V3(:,ind3(K+1:end));
RA3 = V3(:,ind3(1:K))/diag(d3(ind3(1:K)))*V3(:,ind3(1:K))';
% L = size(RS,1);
for i = 1:length(theta1)
    a = exp(-1i*pi*(0:M-1).'*sin(theta1(i)*pi/180));
    P1(i) = 1/(a'*(UN1*UN1')*a);
    P2(i) = 1/(a'*(UN2*UN2')*a);
    P3(i) = (a'*RA2*a)/(a'*(UN2*UN2')*a); 
    a = exp(-1i*pi*(0:Kelm-1).'*sin(theta1(i)*pi/180));
    P4(i) = 1/(a'*(UN3*UN3')*a);
    P5(i) = (a'*RA3*a)/(a'*(UN3*UN3')*a); 
end

figure
plot(theta1,abs(P1),theta1,abs(P2),theta1,abs(P3),theta1,abs(P4),theta1,abs(P5))

inds1 = find_peak(P1,K);
inds2 = find_peak(P2,K);
inds3 = find_peak(P3,K);
inds4 = find_peak(P4,K);
inds5 = find_peak(P5,K);



the1 = theta1(inds1)
the2 = theta1(inds2)
the3 = theta1(inds3)
the4 = theta1(inds4)
the5 = theta1(inds5)

for k = 1:K
    a = exp(-1i*pi*(0:M-1).'*sin(the3(k)/180*pi));
    p(k) = abs(1/(a'*RA2*a));
    a = exp(-1i*pi*(0:Kelm-1).'*sin(the3(k)/180*pi));
    p1(k) = abs(1/(a'*RA3*a));
end