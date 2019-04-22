clear;clc;close all;

load('data\RxSignal_plane_1207.mat');
% load('data\selMat6464.mat')
SNR = 0;
[phi,psi] = sort2D(phi,psi,1);
u = phi*N0*d0;
v = psi*lambda/2*R0/d;
x = P(1) + v - u*u0(1);
y =  P(2) + (- (x-P(1))*u0(1)-u)/u0(2);
rcs = 12*ones(I,1);
figure
scatter(u,v,rcs,'k','filled')
grid on
xlabel('$u/{\mathrm{m}}$','Interpreter','latex')
ylabel('$v/{\mathrm{m}}$','Interpreter','latex')
grid on
% colormap('hot')
figure
scatter(x,y,rcs,'k','filled')
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
grid on
A1 = exp(1i*2*pi*(0:N0-1).'*phi.');
% A1(:,1)'*A1(:,1)
A2 = exp(1i*2*pi*(0:Nv-1).'*psi.');
% rank(A2)
A = kr(A2,A1);
% A*A'
S = 12*exp(1i*2*pi*rand(I,1)*(0:K-1));
% P = S*S'/K;
% diag(P)
X0 =  A*S;
X = X0;
% X = awgn(X0,SNR,'measured');
[U,~,~] = svds(X,I);
J1 = kron(eye(Nr),kron(eye(Nt),[eye(N0-1) zeros(N0-1,1)]));
J2 = kron(eye(Nr),kron(eye(Nt),[zeros(N0-1,1) eye(N0-1)]));
[T,Phi1] = eig((J1*U)\(J2*U));
f10 = angle(diag(Phi1)).'/2/pi/Nt;
A = U*T;
Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f10),ones(N0,I)));
A2 = Phi10.*A;
J3 = kron([eye(Nv-1) zeros(Nv-1,1)],eye(N0));
J4 = kron([zeros(Nv-1,1) eye(Nv-1)],eye(N0));
Phi2 = (J3*A2)\(J4*A2);
f2 = angle(diag(Phi2)).'/2/pi;
Phi20 = kr(exp(-1i*2*pi*(0:Nv-1).'*f2),ones(N0,I));
P = kron(eye(Nr),rowPermutateMat(N0,Nt));
A1 = P*Phi20.*A;
J5 = kron(eye(Nr),[eye(N-1) zeros(N-1,1)]);
J6 = kron(eye(Nr),[zeros(N-1,1) eye(N-1)]);
Phi1 = (J5*A1)\(J6*A1);
f1 = angle(diag(Phi1)).'/2/pi/Nt;
Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f1),exp(-1i*2*pi*(0:N0-1).'*f1*Nt)));
A2 = Phi10.*A;
Phi2 = (J3*A2)\(J4*A2);
f2 = angle(diag(Phi2)).'/2/pi;
% tic
% [p1,s1,S1] = U_ESPRIT2D(X,N0,Nv,I,K1,K2,K3,K4);
% toc
% % tic
% % [p2,s2] = Unitary_ESPRIT_2D1115(X,N0,Nv,I);
% % toc
% 
% ms = sort_matrix([p1 s1 S1],'ascend',1);
% p1 = ms(:,1);
% s1 = ms(:,2);
% S1 = ms(:,3);
% 
% % ms = sort_matrix([p2 s2],'ascend',1);
% % p2 = ms(:,1);
% % s2 = ms(:,2);
% % % S2 = ms(:,3);
% 
% u_est = p1*N0*d0;
% v_est = s1*lambda/2*R0/d;
% figure
% scatter(u_est,v_est,S1,'k','filled');
% grid on
% x_est = v_est - u_est*u0(1);
% y_est =  (- x_est*u0(1)-u_est)/u0(2);
% figure
% scatter(x_est+P(1),y_est+P(2),S1,'k','filled');
% grid on
% 
% % u_est2 = p2*N0*d0;
% % v_est2 = s2*lambda/2*R0/d;
% % figure
% % scatter(u_est2,v_est2,12);
% % x_est2 = v_est2 - u_est2*u0(1);
% % y_est2 =  (- x_est2*u0(1)-u_est2)/u0(2);
% % figure
% % scatter(x_est2,y_est2,S1);
% % tic
% % [p2,s2] = PM2D(X,N0,Nv,I);
% % toc
% % [p2,s2] = sort2D(p2,s2,1);
% % f2 = -0.2:0.001:0.2;
% % tic
% % [p3,s3] = RD_MUSIC(X,f2,N0,Nv,I);
% % toc
% % [p3,s3] = sort2D(p3.',s3.',1);
