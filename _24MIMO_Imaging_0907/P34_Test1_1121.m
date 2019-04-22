clear;clc;close all;


load('data\RxSignal_plane_1223.mat');
ms = sort_matrix([phi,psi],'ascend',1);
phi = ms(:,1);
psi = ms(:,2);
figure(1)
scatter(u,v);
hold on
% Nt = 8;
% Nv = Nt*Nr;
% I = 10;
% I = 8;
% phi = 0.01+(0:I-1)*0.005;
% psi = sin((-30+5*(0:I-1))*pi/180);

A1 = exp(1i*2*pi*(0:N0-1).'*phi.'*Nt);
A2 = exp(1i*2*pi*(0:Nt-1).'*(phi+psi).');
A3 = exp(1i*2*pi*(0:Nr-1).'*psi.'*Nt);
A = kr(A3,kr(A2,A1));
% tmp = kr(A1,A2.*exp(-1i*2*pi*(0:Nt-1).'*(psi(1:I)).'));
% Phi20 = kr(exp(-1i*2*pi*(0:Nr-1).'*psi.'*Nt),kr(exp(-1i*2*pi*(0:Nt-1).'*psi.'),ones(N0,I)));
% tmp1 = Phi20.*A;
% P = kron(eye(Nr),rowPermutateMat(N0,Nt));
% A1 = P*(Phi20.*A);
% J5 = kron(eye(Nr),[eye(N-1) zeros(N-1,1)]);
% J6 = kron(eye(Nr),[zeros(N-1,1) eye(N-1)]);
% Phi1 = (J5*A1)\(J6*A1); 
% % f1 = angle(diag(Phi1)).'/2/pi;
% % max(max(kr(A1,A2) - rowPermutateMat(N0,Nt)*kr(A2,A1)))
% % K = 1;
S = 12*exp(1i*rand(I,1)*(0:K-1));
% % P = S*S'/K;
Y0 = A*S;
SNR = -5;
Y = awgn(Y0,SNR,'measured');

NFFT = 64;
Y1 = reshape(Y(:,1),[N0,Nv]);
for nv = 1:Nv
    X11(:,nv) = fftshift(fft(Y1(:,nv),NFFT))/Nv;
end
for nf = 1:NFFT
    X111(nf,:) = fftshift(fft(X11(nf,:),NFFT))/N0;
end
u_ = (-(NFFT/2-1):NFFT/2)*d0*N0/NFFT;
v_ = (-(NFFT/2-1):NFFT/2)/Nv*lambda*R0/2/d*Nv/NFFT;

figure
imagesc(v_,u_,abs(X111))
axis xy
xlabel('$v/{\mathrm{m}}$','Interpreter','latex')
ylabel('$u/{\mathrm{m}}$','Interpreter','latex')
colormap('gray')

[m,n,S_1] = find_peak_2D(abs(X111),I);
u_1 = u_(m);
v_1 = v_(n);
m_s = sort_matrix([u_1.'+v_1.',u_1.',v_1.',S_1.'],'ascend',1);
u_1 = m_s(:,2);
v_1 = m_s(:,3);
S_1 = m_s(:,4);
x_1 = v_1 - u_1*u0(1);
y_1 = (- x_1*u0(1)-u_1)/u0(2);

figure
scatter(x_1,y_1,mag2db(S_1),'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
[f1,f2,S1] = Uni_Succ_ESPRIT(Y,N0,Nt,Nr,I);
m_s = sort_matrix([f1,f1+f2,f2,S1],'ascend',1);
f1 = m_s(:,1);
f2 = m_s(:,3);
S1 = m_s(:,4);
u1 = f1*N*d0;
v1 = f2*lambda*R0/d/2;
x1 = v1 - u1*u0(1);
y1 =  (- x1*u0(1)-u1)/u0(2);
figure
% hold on
scatter(x1,y1,mag2db(S1),'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')

figure
scatter(x,y,mag2db(targetRCS),'k','filled')
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
figure
scatter(v,u,mag2db(targetRCS),'k','filled')
grid on
xlabel('$v/{\mathrm{m}}$','Interpreter','latex')
ylabel('$u/{\mathrm{m}}$','Interpreter','latex')
% [f11,f12] = ESPRIT_2D(Y,N0,Nv,I);
% [f11,f12] = sort2D(f11/Nt,f12,1);
% % [f1,f2,f3] = MDF_3D(reshape(Y,[N0,Nt,Nr]),I);
% % ms = sort_matrix([f1/Nt,f2,f3],'ascend',1);
% % f1 = ms(:,1);
% % f2 = ms(:,2);
% % f3 = ms(:,3);
% [U,~,~] = svds(Y,I);
% % Te = kron(eye(Nr),rowPermutateMat(N0,Nt));
% J1 = kron(eye(Nv),[eye(N0-1) zeros(N0-1,1)]);
% J2 = kron(eye(Nv),[zeros(N0-1,1) eye(N0-1)]);
% [T,Phi1] = eig((J1*(U))\(J2*(U)));
% f10 = phase(diag(Phi1))/2/pi/Nt;
% [f10,inds] = sort(f10);
% T = T(:,inds);
% A_ = U*T;
% Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f10.'),exp(-1i*2*pi*(0:N0-1).'*f10.'*Nt)));
% A2_ = Phi10.*A_;
% J3 = kron([eye(Nv-1) zeros(Nv-1,1)],eye(N0));
% J4 = kron([zeros(Nv-1,1) eye(Nv-1)],eye(N0));
% Phi2 = (J3*A2_)\(J4*A2_);
% f2 = phase(diag(Phi2))/2/pi;
% f1 = f10;
% % Phi20 = kr(exp(-1i*2*pi*(0:Nr-1).'*f2.'*Nt),kr(exp(-1i*2*pi*(0:Nt-1).'*f2.'),ones(N0,I)));
% % P = kron(eye(Nr),rowPermutateMat(N0,Nt));
% % A1_ = P*(Phi20.*A_);
% % J5 = kron(eye(Nr),[eye(N-1) zeros(N-1,1)]);
% % J6 = kron(eye(Nr),[zeros(N-1,1) eye(N-1)]);
% % Phi1 = (J5*A1_)\(J6*A1_); 
% % f1 = phase(diag(Phi1))/2/pi;
% % Phi10 = kr(ones(Nr,I),kr(exp(-1i*2*pi*(0:Nt-1).'*f1.'),exp(-1i*2*pi*(0:N0-1).'*f1.'*Nt)));
% % A2_ = Phi10.*A_;
% % Phi2 = (J3*A2_)\(J4*A2_);
% % f2 = phase(diag(Phi2))/2/pi;
% u1 = f1*N*d0;
% v1 = f2*lambda*R0/d/2;
% figure(1)
% scatter(u1,v1,'*')
% hold on
% % 
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