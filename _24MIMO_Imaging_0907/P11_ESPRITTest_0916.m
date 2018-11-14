clear;clc;close all;
P = 8;
Q = 6;
d = 2;
SNR = 10;
G = exp(1i*pi*0.1*(0:P-1)).'*exp(1i*pi*0.2*(0:Q-1)) + exp(1i*pi*0.3*(0:P-1)).'*exp(1i*pi*0.4*(0:Q-1));
G = awgn(G,SNR,'measured');
PI_P = IEM(P);
PI_Q = IEM(Q);
PI_2Q = IEM(2*Q);
J1 = [eye(P-1),zeros(P-1,1)];
J2 = rot90(J1,2);
% J2 - IEM(P-1)*J1*IEM(P)
G1 = G(1:P/2,:);
G2 = G(P/2+1:end,:);
GT = [G PI_P*conj(G)*PI_Q];
% IEM(P/2)*conj(G1)*PI_Q - GT(P/2+1:end,Q+1:end)
Q_P = UniMat(P);
Q_2Q = UniMat(2*Q);
GR = real(Q_P'*GT*Q_2Q);
[U,S,V] = svd(GR);
Us = U(:,1:d);
Phi = pinv(Us(1:d,:))*Us(P-d+1:end,:);
[V1 D] = eig(Phi)
theta = atan(D)/pi

% % real(G1+IEM(P/2)*G2) - GR(1:P/2,1:Q)
% % -imag(G1+IEM(P/2)*G2) - GR(1:P/2,1+Q:end)
% % imag(G1-IEM(P/2)*(G2)) - GR(P/2+1:end,1:Q)
% % % real(G1-IEM(P/2)*G2) - GR(P/2+1:end,1+Q:end)
% clear;clc;close all;
% M = 8;
% N = 6;
% aM = exp(1i*pi*0.1*(-M/2:M/2));
% aM = aM(:);
% aN = exp(1i*pi*0.2*(-N/2:N/2));
% aN = aN(:);
% A = aM*aN.';
% QM = UniMat(M+1);
% QN = UniMat(N+1);
% D = QM'*A*QN;

