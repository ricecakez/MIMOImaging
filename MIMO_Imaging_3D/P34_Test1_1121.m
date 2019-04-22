clear;clc;close all;


load('data\RxSignal_plane_1207.mat');
% I = 10;
f = phi+psi;
ms = sort_matrix([phi.',f.',psi.'],'ascend',1);
f = ms(1:I,2);
phi = ms(1:I,1);
psi = ms(1:I,3);
u = phi*N*d0;
v = psi*lambda/2*R0/d;
figure
scatter(u,v,12,'o')
hold on


A1 = exp(1i*2*pi*(0:N0-1).'*phi(1:I).'*Nt);
A2 = exp(1i*2*pi*(0:Nt-1).'*(phi(1:I)+psi(1:I)).');
A3 = exp(1i*2*pi*(0:Nr-1).'*psi(1:I).'*Nt);
A = kr(A3,kr(A2,A1));
% tmp = kr(A1,A2.*exp(-1i*2*pi*(0:Nt-1).'*(psi(1:I)).'));
% Phi20 = kr(exp(-1i*2*pi*(0:Nr-1).'*psi.'*Nt),kr(exp(-1i*2*pi*(0:Nt-1).'*psi.'),ones(N0,I)));
% tmp1 = Phi20.*A;
% P = kron(eye(Nr),rowPermutateMat(N0,Nt));
% A1 = P*(Phi20.*A);
% J5 = kron(eye(Nr),[eye(N-1) zeros(N-1,1)]);
% J6 = kron(eye(Nr),[zeros(N-1,1) eye(N-1)]);
% Phi1 = (J5*A1)\(J6*A1); 
% f1 = angle(diag(Phi1)).'/2/pi;
% max(max(kr(A1,A2) - rowPermutateMat(N0,Nt)*kr(A2,A1)))
S = 12*exp(1i*rand(I,1)*(0:K-1));
P = S*S'/K;
Y = A*S;
SNR= 100;
Y = awgn(Y,SNR,'measured');
[U,~,~] = svds(Y,I);
J1 = kron(eye(Nr),kron(eye(Nt),[eye(N0-1) zeros(N0-1,1)]));
J2 = kron(eye(Nr),kron(eye(Nt),[zeros(N0-1,1) eye(N0-1)]));
[T,Phi1] = eig((J1*U)\(J2*U));
f10 = angle(diag(Phi1)).'/2/pi/Nt;
[f10,inds] = sort(f10);
T = T(:,inds);
A_ = U*T;
Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f10),exp(-1i*2*pi*(0:N0-1).'*f10*Nt)));
A2_ = Phi10.*A_;
J3 = kron([eye(Nv-1) zeros(Nv-1,1)],eye(N0));
J4 = kron([zeros(Nv-1,1) eye(Nv-1)],eye(N0));
Phi2 = (J3*A2_)\(J4*A2_);
f2 = angle(diag(Phi2)).'/2/pi;
f1 = f10;
% Phi20 = kr(exp(-1i*2*pi*(0:Nr-1).'*f2*Nt),kr(exp(-1i*2*pi*(0:Nt-1).'*f2),ones(N0,I)));
% P = kron(eye(Nr),rowPermutateMat(N0,Nt));
% A1_ = P*(Phi20.*A_);
% J5 = kron(eye(Nr),[eye(N-1) zeros(N-1,1)]);
% J6 = kron(eye(Nr),[zeros(N-1,1) eye(N-1)]);
% Phi1 = (J5*A1_)\(J6*A1_); 
% f1 = angle(diag(Phi1)).'/2/pi;
% Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f1),exp(-1i*2*pi*(0:N0-1).'*f1*Nt)));
% A2_ = Phi10.*A_;
% Phi2 = (J3*A2_)\(J4*A2_);
% f2 = angle(diag(Phi2)).'/2/pi;
u1 = f1*N*d0;
v1 = f2*lambda*R0/d/2;
figure(1)
scatter(u1,v1,'*')

% u11 = f11*N/Nt*d0;
% v11 = f12*lambda*R0/d/2;
% scatter(u11,v11,'s')
% [U,D,~] = svd(Y);
% Es = U(:,1:I)*D(1:I,1:I);
% En = U(:,I+1:end);
% Pen = En*En';
% J1 = kron(eye(Nr*Nt),selectMat(N0,1));
% J2 = kron(eye(Nr*Nt),selectMat(N0,2));
% J3 = kron(kron(eye(Nr),selectMat(Nt,1)),eye(N0));
% J4 = kron(kron(eye(Nr),selectMat(Nt,2)),eye(N0));
% Psi_12 = pinv(J3*Es)*(J4*Es);
% [T,Phi_12] = eig(Psi_12);
% [f12,inds] = sort(angle(diag(Phi_12))/2/pi);
% T = T(:,inds);
% A = Es*T;
% Phi_mu = pinv(J1*A)*(J2*A);
% mu = angle(diag(Phi_mu))/2/pi/Nt;
% nu = f12 - mu;
% dt = 0.01;
% delta_dt = 0.001;
% e1 = [1;zeros(Nr*Nt-1,1)];
% for i = 1:I
%     u = (u_ini(i) - dt): delta_dt:(u_ini(i)+dt);
%     for mm = 1:length(u)
%         au = kron(exp(1i*2*pi*(0:N0-1).'*u(mm)*Nt),eye(Nr*Nt));
%         p(mm) = e1'*inv(au'*Pen*au)*e1;
%     end
%     mu_est(i) = u(find(p==max(p)));
% end