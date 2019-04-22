clear;clc;close all;

load('data\RxSignal_plane_0110.mat');

SNR = -5;
% Echo = RxSig/sqrt(N);
Echo = awgn(RxSig/sqrt(N),SNR,'measured');
Ri = rangeangle(Q);
T_0 = 2*R0/c;
% (Ri - R0)/d0;
% Q(1,:)./Ri/2*dt
for k = 1:1
    t = T_min + (k-1)*T + (tc:1/fs:(T-1/fs));
    t = t(:);
    for nr = 1:Nr
        yc = exp(-1i*2*pi*Fc*t).*Echo((k-1)*N_T+((N_c+1):N_T),nr);
        tmp = fft(yc)./ank(:,k)/sqrt(N);
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
             Rcomp(nv) = rangeangle(P,TantPos(:,nt))+rangeangle(RantPos(:,nr),P);
%             Rcomp(nv) = 2*R0-P(1)/R0*dt*(nv-1)+(P(2)^2+P(3)^2)/R0^3*dt^2*((nt-1)^2+((nr-1)*Nt)^2)/2;
%             tau_c = (2*R0-P(1)/R0*dt*(nv-1)...
%                 +1/2*(P(2)^2+P(3)^2)/R0^3*dt^2*((nt-1)^2+((nr-1)*Nt)^2))/c;%...
            %             +(P(2)^2+P(3)^2)/R0^3*dt^2*(nv-1)^2)
            %             tau_c = 1/c*(P(1)/R0*dt*(nv-1)-1/2*(1/R0-P(1)^2/R0^3)*dt^2*((nt-1)^2+((nr-1)*Nt)^2));
            Phi0 = exp(-1i*2*pi*((nt:Nt:N)-1).'*df*T_min).*exp(1i*2*pi*((nt:Nt:N)-1).'/N*Rcomp(nv)/2/d0);
            Psi0 = exp(1i*2*pi*Rcomp(nv)/lambda);
            y(:,nv,k) = tmp(nt:Nt:N).*Phi0*Psi0;
            Y(:,nv,k) = fftshift(fft(tmp(nt:Nt:N).*Phi0*Psi0,N));
        end
    end
end
 for nv = 1:Nv
nnn(nv) = find(Y(:,nv) == max(Y(:,nv)));
 end
figure
plot(nnn)
Rn = (-N/2:(N/2-1))/N*d0*N0;
dR1 = Rn(2)-Rn(1);
 
(((Q(1,1)-P(1)))/Ri(1)*dt)/lambda
figure
imagesc((0:Nv-1),Rn,mag2db(abs(Y)/max(max(abs(Y)))),[-20,0]);
M = 4*Nv;
rho = (-M/2:(M/2-1))/M;
G = zeros(length(Rn),length(rho));
for n = 1:N
    Pc = kron(ones(Nr,1),exp(-1i*2*pi*(0:Nt-1).'/N*Rn(n)/d0)).';
    for nv = 1:Nv
        nt = mod(nv-1,Nt);
        nr = floor((nv-1)/Nt);
        Pc1(nv) = exp(-1i*2*pi*Rn(n)*dt^2*(nt^2+(nr*Nt)^2)/2/R0^2/lambda);
    end
    Z(n,:) = fftshift(fft(Y(n,:),M));
    Z1(n,:) = fftshift(fft(Y(n,:).*Pc.*Pc1,M));
end
figure
imagesc(rho*lambda,Rn,mag2db(abs(Z)/max(max(abs(Z)))),[-20,0])
colormap('hot')
axis xy
figure
imagesc(rho*lambda,Rn,mag2db(abs(Z1)/max(max(abs(Z1)))),[-20,0])
colormap('hot')
axis xy
% ylim([-10,20])
R1 = R0-Ri;
qi = Q - P;
mu1 = (qi(1,:)+R1*u0(1))/R0*dt;
figure
scatter(mu1,R1)

% n = 539;
for m = 1:M
    for n = 1:N
        for nv = 1:Nv
            Rnm(nv) =  rho(m)*lambda*(nv-1)/2 + Rn(n);
            nn(m,nv) = round((Rnm(nv)-Rn(1))/dR1)+1;
            if nn(m,nv) < 1
                nn(m,nv) = N+nn(m,nv);
            end
            if nn(m,nv) > N
                nn(m,nv) = nn(m,nv)-N;
            end
            nt = mod(nv-1,Nt);
            nr = (nv-nt-1)/Nt;
%             if nt == 0
%                 nt = Nt;
%             end
            G(n,m) = G(n,m) + Y(nn(m,nv),nv)*exp(-1i*2*pi*nt/2/N*Rnm(nv)/d0)...
                *exp(-1i*2*pi*rho(m)*(nv-1))*exp(-1i*2*pi*Rn(n)/2/R0^2/lambda*dt^2*(nt^2+(nr*Nt)^2));     
        end
    end    
end
% [n1,m1] = find(G==max(max(G)))
% [n2,m2] = find(Z==max(max(Z)))
% figure
% plot(abs(G(n1,:)))
% hold on
% plot(abs(Z(n2,:)))
figure
imagesc(rho,Rn,mag2db(abs(G)/max(max(abs(G)))),[-20,0])
colormap('hot')
axis xy
[m,n,S_1] = find_peak_2D(abs(G),I);
R_1 = Rn(m);
rho_1 = rho(n);
figure(5)
hold on
scatter(rho_1*lambda,R_1,'*')
[m,n,S_1] = find_peak_2D(abs(Z1),I);
R_1 = Rn(m);
rho_1 = rho(n);
figure(5)
hold on
scatter(rho_1*lambda,R_1,'+')
% figure
% imagesc(v_,u_,mag2db(abs(Y)/max(max(abs(Y)))),[-30 0])
% figure
% imagesc(v_,u_,mag2db(abs(Y1)/max(max(abs(Y1)))),[-30,0])
% figure
% imagesc(v_,u_,mag2db(abs(Y2)/max(max(abs(Y2)))),[-30,0])

