clear;clc;close all;
load('data\RxSignal_plane_1207.mat');
phi = phi.';
psi1 = psi1.';
psi2 = psi2.';
f = phi+psi1;
ms = sort_matrix([phi,f,psi1,psi2,a.'],'ascend',1);
f = ms(:,2);
phi = ms(:,1);
psi1 = ms(:,3);
% f - phi - psi1
psi2 = ms(:,4);
x = ms(:,5);
y = ms(:,6);
z = ms(:,7);
% Nt = 16;
% N0 = 16;
% Nr = 16;
% I = 15;
% phi = -0.04 + 0.005*(0:I-1);
% psi1 = -0.45 + 0.06*(0:I-1);
% psi2 = -0.38 + 0.04*(0:I-1);
% f = phi+psi1;
u = phi*N*d0;
v1 = psi1*lambda/2*R0/dt;
v2 = psi2*lambda/2*R0/dr;
figure
scatter3(a(1,:),a(2,:),a(3,:),targetRCS*10,'filled')
hold on
xlabel('x')
ylabel('y')
zlabel('z')
A1 = exp(1i*2*pi*(0:N0-1).'*phi.'*Nt);
% phi*Nt
A2 = exp(1i*2*pi*(0:(Nt-1)).'*f.');
A3 = exp(1i*2*pi*(0:(Nr-1)).'*psi2.');
% A = kr(A3,kr(A2,A1));
% C = diag(targetRCS);

G1 = [ones(N0,1) (0:N0-1).'];
G2 = [ones(Nt,1) (0:Nt-1).'];
G3 = [ones(Nr,1) (0:Nr-1).'];
% pinv(G3)*phase(A3(:,1))/2/pi
X0 = krb(A2,A1)*A3.';
SNR = 5;
X = awgn(X0,SNR,'measured');
X = tensor(reshape(X,[N0,Nt,Nr]));
% A = double(tenmat(X,1))/(krb(A3,A2)).';
% B = double(tenmat(X,2))/(krb(A3,A1)).';
% C = double(tenmat(X,3))/(krb(A2,A1)).';
% X = full(ktensor({A1,A2,A3}));
[p1,p2,p3] = MDF_3D(X,I);
% X1 = UniTrans(X,1);
% [A,B,C] = TALS(X1,I);
% M = size(A,1);
% N = size(B,1);
% P = size(C,1);
% J12 = [zeros(M-1,1) eye(M-1)];
% tmp = UniMat(M-1)'*J12*UniMat(M);
% K11 = real(tmp);
% K12 = imag(tmp);
% J22 = [zeros(N-1,1) eye(N-1)];
% tmp = UniMat(N-1)'*J22*UniMat(N);
% K21 = real(tmp);
% K22 = imag(tmp);
% J32 = [zeros(P-1,1) eye(P-1)];
% tmp = UniMat(P-1)'*J32*UniMat(P);
% K31 = real(tmp);
% K32 = imag(tmp);
% phi_a = (K11*A)\(K12*A);
% phi_b = (K21*B)\(K22*B);
% phi_c = (K31*C)\(K32*C);
% p1 = atan(diag(phi_a))/pi;
% p2 = atan(diag(phi_b))/pi;
% p3 = atan(diag(phi_c))/pi;
% for i = 1:I
% a1 = phase(f1(:,i))/2/pi;
% tmp = pinv(G1)*a1;
% p1(i) = tmp(2);
% a2 = phase(f2(:,i))/2/pi;
% tmp = pinv(G1)*a2;
% p2(i) = tmp(2);
% a3 = phase(f3(:,i))/2/pi;
% tmp = pinv(G1)*a3;
% p3(i) = tmp(2);
% end
ms = sort_matrix([p1/Nt,p2,p3],'ascend',1);
phi_est = ms(:,1);
psi1_est = ms(:,2) - phi_est;
psi2_est = ms(1:I,3);
u_est = phi_est*N*d0;
v1_est = psi1_est*lambda/2*R0/dt;
v2_est = psi2_est*lambda/2*R0/dr;
x_est = v1_est - u_est*u0(1);
y_est = v2_est - u_est*u0(2);
z_est = (-u_est - x_est*u0(1) - y_est*u0(2))/u0(3);
% figure
scatter3(x_est,y_est,z_est,'filled')