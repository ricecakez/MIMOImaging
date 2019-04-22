clear;clc;close all;
load('data\targetPara_0422.mat');
% open('result\fig7_0422.fig')
hold on
c = 3e8;
% SNR = 5;
[f0,lambda,Bw,fs,N,K,df,tb,tc,ts,PRF,Nt,Nr,dt,dr] = MIMORadarPara1;
deltaR = c/2/Bw;
L_T = fs*ts;
L_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
Nv = Nt*Nr;
% [rho,SP,V,P] = TargetPara;
I = size(SP,1);
Q = SP.';
% P = P.';
q = Q - P;
% figure
% scatter(q(1,:),q(2,:),20*log10(rho.'),'filled','k')
% grid on
% xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
% ylabel('$y/{\mathrm{m}}$','Interpreter','latex')

T = [dt*(0:Nt-1)-dt*(Nt-1)/2;zeros(2,Nt)];
R = [dr*(0:Nr-1)-dr*(Nr-1)/2;zeros(2,Nr)];
OP = rangeangle(P);
u0 = P/OP;
ui = q.'*u0;
vi = (q-u0*ui.').'*[1;0;0];
f1i = -ui.'/N/deltaR;
f2i = vi.'*dt/lambda/OP;
% dqi = rangeangle(q);
% M_sort = sort_matrix([f1i.'+f2i.';f1i.';f2i.'],'ascend',2);
% f1i = M_sort(2,:);
% f2i = M_sort(3,:);
% ui = -f1i*N*deltaR;
% vi = f2i/dt*lambda*OP;
% xi = vi+ui*u0(1);
% yi = (ui-xi*u0(1))/u0(2);
% figure
% scatter(ui,vi,20*log10(rho.'),'filled','k')
% grid on
% xlabel('$u/{\mathrm{m}}$','Interpreter','latex')
% ylabel('$v/{\mathrm{m}}$','Interpreter','latex')
C1 = exp(1i*2*pi*(0:N0-1).'*Nt*f1i(1:I));
C2 = exp(1i*2*pi*(0:Nv-1).'*f2i(1:I));
qt = ui(1:I).'/OP/OP/lambda*dt^2;
Dr = exp(1i*pi*((0:Nr-1).').^2*qt*Nt^2);
Dt = exp(1i*pi*((0:Nt-1).').^2*qt)...
    .*exp(1i*2*pi*(0:Nt-1).'*f1i);
D = krb(Dr,Dt);
B = diag(rho)*exp(1i*randn(I,K));
B1 = diag(exp(-1i*4*pi*ui/lambda))*B;
C = krb((D.*C2),C1);
% J1 = kron(eye(Nv),[eye(N0-1),zeros(N0-1,1)]);%*kron(eye(Nr),permuMat(N0,Nt));
% J2 = kron(eye(Nv),[zeros(N0-1,1),eye(N0-1)]);%*kron(eye(Nr),permuMat(N0,Nt));
% J3 = kron([eye(Nv-1),zeros(Nv-1,1)],eye(N0));
% J4 = kron([zeros(Nv-1,1),eye(Nv-1)],eye(N0));

m = N0 - 1;
for nv = 1:Nv
    J11((nv-1)*m+(1:m),1) = (nv-1)*N0+(1:m);
    J21((nv-1)*m+(1:m),1) = (nv-1)*N0+(2:N0);
end

m = Nv - 1;
J31 = 1:N0*m;
J41 = (N0+1):N0*Nv;

Y0 = C*B1(1:I,:);
Py = (mean(mean(abs(Y0).^2)));
NFFT = 512;
OO = 100;
% snr_inds = 1;
SNR_db = -5;
% SNR_db = [-20:2:20];
MM = 20;
NN = 150;
eps1 = 0.5e-3;
eps2 = 1e-3;
% for snr_inds = 1:length(SNR_db)
    % SNR_db = -5;
    SNR = 10^(SNR_db/10);  %10db
%     rmse_x_FFT(snr_inds) = 0;
%     rmse_y_FFT(snr_inds) = 0;
%     rmse_x_ce(snr_inds) = 0;
%     rmse_y_ce(snr_inds) = 0;
%     rmse_x_me(snr_inds) = 0;
%     rmse_y_me(snr_inds) = 0;
    %     rmse_x_fe(snr_inds) = 0;
    %     rmse_y_fe(snr_inds) = 0;
    
%     for oo = 1:OO
        wn = (randn(N0*Nv,K) + 1i*randn(N0*Nv,K))*sqrt(Py/SNR/2);
        YY = Y0+wn;%awgn(Y0,SNR_db,'measured');
%         Y1 = reshape(YY(:,1),[N0,Nv]);
%         for n0 = 1:N0
%             F1(n0,:) = fftshift(fft(Y1(n0,:),NFFT))/Nv;
%         end
%         for nf = 1:NFFT
%             F2(:,nf) = fftshift(fft(F1(:,nf),NFFT))/N0;
%         end
%         uf = -linspace(-0.5,0.5,NFFT)*N0*deltaR;
%         vf = linspace(-0.5,0.5,NFFT)*lambda*OP/dt;
%         figure
%         imagesc(uf,vf,20*log10(abs(F2.')))
%         axis xy;
%         xlabel('$u/{\mathrm{m}}$','Interpreter','latex')
%         ylabel('$v/{\mathrm{m}}$','Interpreter','latex')
%         colormap('jet')
%         [mm,nn,Smax] = find_peak_2D(20*log10(abs(F2)),I);
%         uf_est = uf(mm);
%         vf_est = vf(nn);
%         %         [uf_est,vf_est] = sort2D(uf_est,vf_est,2);
%         xf_est = vf_est+uf_est*u0(1);
%         yf_est = (uf_est-xf_est*u0(1))/u0(2);
%         [xf_est1,yf_est1] = sort_test(xf_est,yf_est,q(1,:),q(2,:));
%         rmse_x_FFT(snr_inds) = rmse_x_FFT(snr_inds)+sum((xf_est1-q(1,:)).^2);
%         rmse_y_FFT(snr_inds) = rmse_y_FFT(snr_inds)+sum((yf_est1-q(2,:)).^2);
        % figure
        % scatter(xf_est,yf_est,Smax,'*','k')
        % hold on
        % grid on
        % xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
        % ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
        % hold on
        % scatter(q(1,:),q(2,:),20*log10(rho.'),'r')
        % legend('Reconstructed scatterers using FFT','Original scatterers')
        nnv = Nv;
        [Es] = Uni_ss_est(YY,I);
        Psi = (Es(J11,:))\(Es(J21,:));
        [T1,D1] = eig(Psi);
        f1_est0 = angle(diag(D1).')/2/pi/Nt;
        A = Es*T1;
        
        %% Conventional ESPRIT
        Phi1 = (A(J31,:))\(A(J41,:));
        f2_est0 = angle(diag(Phi1)).'/2/pi;
        ui_est0 = -f1_est0*N*deltaR;
        vi_est0 = f2_est0*lambda*OP/dt;
        xc_est = vi_est0+ui_est0*u0(1);
        yc_est = (ui_est0-xc_est*u0(1))/u0(2);
        
        %÷ÿ–¬≈≈–Ú
        [xc_est1,yc_est1,indsc] = sort_test(xc_est,yc_est,q(1,:),q(2,:));
        ui_est = ([xc_est1.' yc_est1.' zeros(I,1)]*u0).';
        vi_est = (([xc_est1.' yc_est1.' zeros(I,1)].'-u0*ui_est).'*[1;0;0]).';
        f1_est = -ui_est/N/deltaR;
        f2_est = vi_est*dt/lambda/OP;
        C1_est = exp(1i*2*pi*(0:N0-1).'*f1_est*Nt);
        C2_est = exp(1i*2*pi*(0:Nv-1).'*f2_est);
        C_est = kr(C2_est,C1_est);
        S_est = 20*log10(mean(abs((C_est'*C_est)\C_est'*YY),2)).';
%         inds = find(S_est < 0);
%         S_est(inds) = 1;
        Dr_est = exp(1i*pi*((0:Nr-1).').^2*ui_est*dr^2/OP^2/lambda);
        Dt_est = exp(1i*pi*((0:Nt-1).').^2*ui_est*dt^2/OP^2/lambda)...
            .*exp(1i*2*pi*(0:Nt-1).'*f1_est);
        D_est = krb(Dr_est,Dt_est);
%         rmse_x_ce(snr_inds) = rmse_x_ce(snr_inds)+sum((xc_est1-q(1,:)).^2);
%         rmse_y_ce(snr_inds) = rmse_y_ce(snr_inds)+sum((yc_est1-q(2,:)).^2);
        %         Comp = kron(conj(D_est),ones(N0,1));
        %         A1 = Comp.*A;
        %         Phi = (A1(J31,:))\(A1(J41,:));
        %         f2_est = angle(diag(Phi)).'/2/pi;
        %         M_sort = sort_matrix([f1_est+f2_est;f1_est;f2_est],'ascend',2);
        %         f1_est = M_sort(2,:);
        %         f2_est = M_sort(3,:);
        
        
        Rn = eye(size(Es,1)) - Es*Es';
        for i = 1:I
            for nn = 1:NN
                f2n(i,nn) = f2_est(i) +(nn-1- NN/2)*eps2;
                c2n = exp(1i*2*pi*(0:Nv-1).'*f2n(i,nn));
                cn = kron((D_est(:,i).*c2n),C1_est(:,i)); 
                P_MUSIC(nn) = 1/(cn'*Rn*cn);
            end
            n1 = find(P_MUSIC==max(P_MUSIC));
            clear P_MUSIC
            f2_est1(i) = f2n(i,n1);
            c2i = exp(1i*2*pi*(0:Nv-1).'*f2_est1(i));
            for mm = 1:MM
                f1m(i,mm) = f1_est(i) +(mm-1- MM/2)*eps1;
                qim(i,mm) = -f1m(i,mm)*N*deltaR*dt^2/OP^2/lambda;
                c1m = exp(1i*2*pi*(0:N0-1).'*f1m(i,mm)*Nt);
                drm = exp(1i*pi*((0:Nr-1).').^2*qim(i,mm)*Nt^2);
                dtm = exp(1i*pi*((0:Nt-1).').^2*qim(i,mm))...
                    .*exp(1i*2*pi*(0:Nt-1).'*f1m(i,mm));
                dm = kron(drm,dtm);
                cm = kron((dm.*c2i),c1m);
                P_MUSIC(mm) = 1/(cm'*Rn*cm);
            end
            m1 = find(P_MUSIC==max(P_MUSIC));
            f1_est1(i) = f1m(i,m1);
%             c1_est1 = exp(1i*2*pi*(0:N0-1).'*f1_est1(i)*Nt);
%             qi_est1(i) = -f1_est1(i)*N*deltaR*dt^2/OP^2/lambda;
%             dr_est1 = exp(1i*pi*((0:Nr-1).').^2*qi_est1(i)*Nt^2);
%             dt_est1 = exp(1i*pi*((0:Nt-1).').^2*qi_est1(i))...
%                 .*exp(1i*2*pi*(0:Nt-1).'*f1_est1(i));
%                 d_est1 = kron(dr_est1,dt_est1);
            clear P_MUSIC
            
           
        end
        %         M_sort = sort_matrix([f1_est1+f2_est1;f1_est1;f2_est1],'ascend',2);
        %         f1_est1 = M_sort(2,:);
        %         f2_est1 = M_sort(3,:);
        
        % [f1,f2] = Uni_ESPRIT_2D(Y2,N0,Nv,I);
        
        %         vi_est = f2_est*lambda*OP/dt;
        ui_est1 = -f1_est1*N*deltaR;
        vi_est1 = f2_est1*lambda*OP/dt;
        
        
        % figure
        % scatter(xi_est0,yi_est0,S_est0,'*','k')
        % hold on
        % grid on
        % xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
        % ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
        % % hold on
        % scatter(q(1,:),q(2,:),20*log10(rho.'),'r')
        % legend('Reconstructed scatterers using conventional 2D-ESPRIT','Original scatterers')
        
        %         %% Modified ESPRIT
        %         xi_est = vi_est+ui_est*u0(1);
        %         yi_est = (ui_est-xi_est*u0(1))/u0(2);
        %         C1_est = exp(1i*2*pi*(0:N0-1).'*f1_est*Nt);
        %         C2_est = exp(1i*2*pi*(0:Nv-1).'*f2_est);
        %         qt_est = ui_est(1:I)/OP/OP/lambda*dt^2;
        %         Dr_est = exp(1i*pi*((0:Nr-1).').^2*qt*Nt^2);
        %         Dt_est = exp(1i*pi*((0:Nt-1).').^2*qt)...
        %             .*exp(1i*2*pi*(0:Nt-1).'*f1_est);
        %         D_est = krb(Dr_est,Dt_est);
        %         C_est = kr((D_est.*C2_est),C1_est);
        %         S_est = 20*log10(mean(abs((C_est'*C_est)\C_est'*YY),2)).';
%                 rmse_x_me(snr_inds) = rmse_x_me(snr_inds)+sum((xi_est-xi).^2);
%                 rmse_y_me(snr_inds) = rmse_y_me(snr_inds)+sum((yi_est-yi).^2);
%         % figure
        % scatter(xi_est,yi_est,S_est,'*','k')
        % hold on
        % grid on
        % xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
        % ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
        % % hold on
        % scatter(q(1,:),q(2,:),20*log10(rho.'),'r')
        % legend('Coarse estimation after modified 2D-ESPRIT','Original scatterers')
        
        %% Modified ESPRIT-MUSIC
        xm_est = vi_est1+ui_est1*u0(1);
        ym_est = (ui_est1-xm_est*u0(1))/u0(2);
        [xm_est1,ym_est1,indsm] = sort_test(xm_est,ym_est,q(1,:),q(2,:));
        ui_est11 = ([xm_est1.' ym_est1.' zeros(I,1)]*u0).';
        vi_est11 = (([xm_est1.' ym_est1.' zeros(I,1)].'-u0*ui_est11).'*[1;0;0]).';
        f1_est11 = -ui_est11/N/deltaR;
        f2_est11 = vi_est11*dt/lambda/OP;
        C1_est1 = exp(1i*2*pi*(0:N0-1).'*f1_est11*Nt);
        C2_est1 = exp(1i*2*pi*(0:Nv-1).'*f2_est11);
        qt_est1 = ui_est11(1:I)/OP/OP/lambda*dt^2;
        Dr_est1 = exp(1i*pi*((0:Nr-1).').^2*qt_est1*Nt^2);
        Dt_est1 = exp(1i*pi*((0:Nt-1).').^2*qt_est1)...
            .*exp(1i*2*pi*(0:Nt-1).'*f1_est11);
        D_est1 = krb(Dr_est1,Dt_est1);
        C_est1 = kr((D_est1.*C2_est1),C1_est1);
        S_est1 = 20*log10(mean(abs((C_est1'*C_est1)\C_est1'*YY),2)).';
        inds = find(S_est1<0);
        S_est1(inds) = 1;
%         rmse_x_me(snr_inds) = rmse_x_me(snr_inds)+sum((xm_est1-q(1,:)).^2);
%         rmse_y_me(snr_inds) = rmse_y_me(snr_inds)+sum((ym_est1-q(2,:)).^2);
%     end
% end
% rmse_x_FFT1 = sqrt(rmse_x_FFT/OO)/I;
% rmse_y_FFT1 = sqrt(rmse_y_FFT/OO)/I;
% rmse_x_ce1 = sqrt(rmse_x_ce/OO)/I;
% rmse_y_ce1 = sqrt(rmse_y_ce/OO)/I;
% rmse_x_me1 = sqrt(rmse_x_me/OO)/I;
% rmse_y_me1 = sqrt(rmse_y_me/OO)/I;
% figure
% semilogy(SNR_db,rmse_x_FFT1,'-+',SNR_db,rmse_x_ce1,'-o',SNR_db,rmse_x_me1,'-*')
% xlabel('${\mathrm{SNR}}/{\mathrm{dB}}$','Interpreter','latex')
% ylabel('${\mathrm{RMSE}}/{\mathrm{m}}$','Interpreter','latex')
% figure
% semilogy(SNR_db,rmse_x_FFT1,SNR_db,rmse_x_ce1,SNR_db,rmse_x_me1)
% xlabel('$SNR/{\mathrm{dB}}$','Interpreter','latex')
% rmse_x_fe1 = sqrt(rmse_x_me/OO)/I;
% rmse_y_fe1 = sqrt(rmse_y_fe/OO)/I;
figure
scatter(q(1,:),q(2,:),20*log10(rho.'),'k')
hold on
% scatter(xf_est,yf_est,Smax,'+','b')
% scatter(xc_est1,yc_est1,S_est,'x','r')
scatter(xm_est1,ym_est1,S_est1,'+','b')
% 
grid on
xlabel('$x/{\mathrm{m}}$','Interpreter','latex')
ylabel('$y/{\mathrm{m}}$','Interpreter','latex')
% hold on

legend('Original scatterers','With $\Delta f = 0.5\mathrm{MHz}$')

% ui_est = -angle(phi_est)/2/pi*N*deltaR;
% figure
% hold on
% scatter(ui,vi)
% scatter(ui_est,vi_est0,'+')
% scatter(ui_est,vi_est,'*')
% scatter(ui_est1,vi_est1,'x')

