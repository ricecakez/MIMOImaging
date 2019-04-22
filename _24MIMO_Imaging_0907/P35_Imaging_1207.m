clear;clc;close all;

load('data\RxSignal_plane_1223.mat');
% Echo = RxSig/sqrt(N);

SNR = -20;
MM = 500;
ms = sort_matrix([phi,psi],'ascend',1);
phi = ms(:,1);
psi = ms(:,2);
u = phi*N*d0;
v = psi*lambda*R0/d/2;
x = v - u*u0(1);
y =  (- x*u0(1)-u)/u0(2);
u_ = (-(N0/2-1):N0/2)*d0;
v_ = (-(Nv/2-1):Nv/2)/Nv*lambda*R0/2/d;

rmse_fft_x = zeros(size(SNR));
rmse_fft_y = zeros(size(SNR));
rmse_succ_uni_x = zeros(size(SNR));
rmse_succ_uni_y = zeros(size(SNR));
rmse_succ_x = zeros(size(SNR));
rmse_succ_y = zeros(size(SNR));
rmse_OFDM_uni_x = zeros(size(SNR));
rmse_OFDM_uni_y = zeros(size(SNR));
rmse_OFDM_x = zeros(size(SNR));
rmse_OFDM_y = zeros(size(SNR));
rmse_2D_uni_x = zeros(size(SNR));
rmse_2D_uni_y = zeros(size(SNR));
rmse_2D_x = zeros(size(SNR));
rmse_2D_y = zeros(size(SNR));

for nn = 1:length(SNR)
    for mm = 1:MM
        
        %         Y = awgn(Y0,SNR(nn),'measured');
        Echo = awgn(RxSig/sqrt(N),SNR(nn),'measured');
        for k = 1:K
            t = T_min + (k-1)*T + (tc:1/fs:(T-1/fs));
            t = t(:);
            for nr = 1:Nr
                yc = exp(-1i*2*pi*Fc*t).*Echo((k-1)*N_T+((N_c+1):N_T),nr);
                tmp = fft(yc)./ank(:,k)/sqrt(N);
                for nt = 1:Nt
                    nv = (nr-1)*Nt+nt;
                    TP = rangeangle(TantPos(:,nt),P);
                    PR = rangeangle(P,RantPos(:,nr));
                    Phi0 = exp(1i*2*pi*((nt:Nt:N)-1).'*df*((TP+PR)/c-T_min));
                    Psi0 = exp(1i*2*pi*(TP+PR)/lambda);
                    Y((nv-1)*N0+(1:N0),k) = tmp(nt:Nt:N).*Phi0*Psi0;
                end
            end
        end
        
        % NFFT = 64;
        % X111 = zeros(NFFT);
        
        % for k = 1:K
        Y1 = reshape(Y(:,1),[N0,Nv]);
        for nv = 1:nv
            X11(:,nv) = fftshift(ifft(Y1(:,nv)));
        end
        for nf = 1:N0
            X111(nf,:) = fftshift(ifft(X11(nf,:)));
        end
        
        
        % imagesc(u_,v_,abs(X111))
        % axis xy
        % xlabel('$u/{\mathrm{m}}$','Interpreter','latex')
        % ylabel('$v/{\mathrm{m}}$','Interpreter','latex')
        % colormap('hot')
        [m,n,S_1] = find_peak_2D(abs(X111),I);
        u_1 = -u_(m);
        v_1 = -v_(n);
        m_s = sort_matrix([u_1.',v_1.',S_1.'],'ascend',1);
        u_1 = m_s(:,1);
        v_1 = m_s(:,2);
        S_1 = m_s(:,3);
        x_1 = v_1 - u_1*u0(1);
        y_1 = (- x_1*u0(1)-u_1)/u0(2);
        
        rmse_fft_x(nn) = rmse_fft_x(nn) + sum((x_1-x).^2);
        rmse_fft_y(nn) = rmse_fft_y(nn) + sum((y_1-y).^2);
        % figure
        % scatter(x_1,y_1,S_1,'s');
        % hold on
        % figure
        % scatter(u_1,v_1,S_1,'s');
        % hold on
        % % [U,~,~] = svds(Y,I);
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
        % ms = sort_matrix([f1,f2],'ascend',1);
        % f1 = ms(:,1);
        % f2 = ms(:,2);
        [f1,f2,S1] = Uni_Succ_ESPRIT(Y,N0,Nt,Nr,I);
        m_s = sort_matrix([f1,f2,S1],'ascend',1);
        f1 = m_s(:,1);
        f2 = m_s(:,2);
        S1 = m_s(:,3);
        u1 = f1*N*d0;
        v1 = f2*lambda*R0/d/2;
        x1 = v1 - u1*u0(1);
        y1 =  (- x1*u0(1)-u1)/u0(2);
        rmse_succ_uni_x(nn) = rmse_succ_uni_x(nn) + sum((x1-x).^2);
        rmse_succ_uni_y(nn) = rmse_succ_uni_y(nn) + sum((y1-y).^2);
        [f12,f22,S12] = Succ_ESPRIT(Y,N0,Nt,Nr,I);
        m_s = sort_matrix([f12,f22,S12],'ascend',1);
        f12 = m_s(:,1);
        f22 = m_s(:,2);
        S12 = m_s(:,3);
        u12 = f12*N*d0;
        v12 = f22*lambda*R0/d/2;
        x12 = v12 - u12*u0(1);
        y12 =  (- x12*u0(1)-u12)/u0(2);
        rmse_succ_x(nn) = rmse_succ_x(nn) + sum((x12-x).^2);
        rmse_succ_y(nn) = rmse_succ_y(nn) + sum((y12-y).^2);
        % figure(3)
        % scatter(u1,v1,'*')
        % hold on
        % scatter(u,v)
        % figure(2)
        % scatter(x1,y1,'*')
        % hold on
        % scatter(x,y)
        %
        [f11,f21,S11] = Uni_ESPRIT_OFDM_2D(Y,N0,Nt,Nr,I);
        m_s = sort_matrix([f11,f21,S11],'ascend',1);
        f11 = m_s(:,1);
        f21 = m_s(:,2);
        S11 = m_s(:,3);
        u11 = f11*N*d0;
        v11 = f21*lambda*R0/d/2;
        x11 = v11 - u11*u0(1);
        y11 =  (- x11*u0(1)-u11)/u0(2);
        rmse_OFDM_uni_x(nn) = rmse_OFDM_uni_x(nn) + sum((x11-x).^2);
        rmse_OFDM_uni_y(nn) = rmse_OFDM_uni_y(nn) + sum((y11-y).^2);
        
        [f112,f122,S112] = ESPRIT_OFDM_2D(Y,N0,Nt,Nr,I);
        m_s = sort_matrix([f112,f122,S112],'ascend',1);
        f112 = m_s(:,1);
        f122 = m_s(:,2);
        S112 = m_s(:,3);
        u112 = f112*N*d0;
        v112 = f122*lambda*R0/d/2;
        x112 = v112 - u112*u0(1);
        y112 =  (- x112*u0(1)-u112)/u0(2);
        rmse_OFDM_x(nn) = rmse_OFDM_x(nn) + sum((x112-x).^2);
        rmse_OFDM_y(nn) = rmse_OFDM_y(nn) + sum((y112-y).^2);
        % f11 = f11/Nt;
        % f12 = f12/Nt;
        [f111,f121,S111] = Uni_ESPRIT_2D(Y,N0,Nv,I);
        m_s = sort_matrix([f111,f121,S111],'ascend',1);
        f111 = m_s(:,1);
        f121 = m_s(:,2);
        S111 = m_s(:,3);
        f111 = f111/Nt;
        
        % ms = sort_matrix([f11/Nt,f12],'ascend',1);
        % f11 = ms(:,1);
        % f12 = ms(:,2);
        
        % figure(3)
        % scatter(u11,v11,'d')
        % hold on
        % scatter(u,v,'*')
        % figure(2)
        % scatter(x11,y11,'d')
        
        u111 = f111*N*d0;
        v111 = f121*lambda*R0/d/2;
        x111 = v111 - u111*u0(1);
        y111 =  (- x111*u0(1)-u111)/u0(2);
        rmse_2D_uni_x(nn) = rmse_2D_uni_x(nn) + sum((x111-x).^2);
        rmse_2D_uni_y(nn) = rmse_2D_uni_y(nn) + sum((y111-y).^2);
        
        [f1112,f1122,S1112] = ESPRIT_2D(Y,N0,Nv,I);
        m_s = sort_matrix([f1112,f1122,S1112],'ascend',1);
        f1112 = m_s(:,1);
        f1122 = m_s(:,2);
        S1112 = m_s(:,3);
        f1112 = f1112/Nt;
        
        u1112 = f1112*N*d0;
        v1112 = f1122*lambda*R0/d/2;
        x1112 = v1112 - u1112*u0(1);
        y1112 =  (- x1112*u0(1)-u1112)/u0(2);
        rmse_2D_x(nn) = rmse_2D_x(nn) + sum((x1112-x).^2);
        rmse_2D_y(nn) = rmse_2D_y(nn) + sum((y1112-y).^2);
    end
    rmse_fft_x(nn) = sqrt(rmse_fft_x(nn)/MM/I);
    rmse_fft_y(nn) = sqrt(rmse_fft_y(nn)/MM/I);
    rmse_succ_uni_x(nn) = sqrt(rmse_succ_uni_x(nn)/MM/I);
    rmse_succ_uni_y(nn) = sqrt(rmse_succ_uni_y(nn)/MM/I);
    rmse_succ_x(nn) = sqrt(rmse_succ_x(nn)/MM/I);
    rmse_succ_y(nn) = sqrt(rmse_succ_y(nn)/MM/I);
    rmse_OFDM_uni_x(nn) = sqrt(rmse_OFDM_uni_x(nn)/MM/I);
    rmse_OFDM_uni_y(nn) = sqrt(rmse_OFDM_uni_y(nn)/MM/I);
    rmse_OFDM_x(nn) = sqrt(rmse_OFDM_x(nn)/MM/I);
    rmse_OFDM_y(nn) = sqrt(rmse_OFDM_y(nn)/MM/I);
    rmse_2D_uni_x(nn) = sqrt(rmse_2D_uni_x(nn)/MM/I);
    rmse_2D_uni_y(nn) = sqrt(rmse_2D_uni_y(nn)/MM/I);
    rmse_2D_x(nn) = sqrt(rmse_2D_x(nn)/MM/I);
    rmse_2D_y(nn) = sqrt(rmse_2D_y(nn)/MM/I);
end
figure
semilogy(SNR,rmse_fft_x,'b-o',SNR,rmse_2D_x,'r--s',SNR,rmse_2D_uni_x,'r-s',...
    SNR,rmse_OFDM_x,'g--*',SNR,rmse_OFDM_uni_x,'g-*',SNR,rmse_succ_x,'k--+',SNR,rmse_succ_uni_x,'k-+')
figure
semilogy(SNR,rmse_fft_y,'b-o',SNR,rmse_2D_y,'r--s',SNR,rmse_2D_uni_y,'r-s',...
    SNR,rmse_OFDM_y,'g--*',SNR,rmse_OFDM_uni_y,'g-*',SNR,rmse_succ_y,'k--+',SNR,rmse_succ_uni_y,'k-+')

% figure(3)
% scatter(u111,v111,'+')
% % hold on
% % scatter(u,v,'*')
% figure(2)
% scatter(x111,y111,'+')
% hold on
% scatter(a(1,:),a(2,:),'*')