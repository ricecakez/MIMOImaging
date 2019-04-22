clear;clc;close all;

load('data\RxSignal_targetPlane_1_0409.mat');

f1i = -ui/N/deltaR;
f2i = vi*dt/lambda/OP;
figure
scatter(ui,vi)
SNR = 1;  %10db
Nv = Nr*Nt; 
% for nr = 1:Nr
Py = sqrt(mean(mean(abs(Y).^2)));
wn = (randn(N*Nr,K) + 1i*randn(N*Nr,K))/sqrt(2)*Py/SNR;
cov(wn);
Yw = Y + wn;
J1 = kron(eye(Nv),[eye(N0-1),zeros(N0-1,1)]);%*kron(eye(Nr),permuMat(N0,Nt));
J2 = kron(eye(Nv),[zeros(N0-1,1),eye(N0-1)]);%*kron(eye(Nr),permuMat(N0,Nt));
J3 = kron([eye(Nv-1),zeros(Nv-1,1)],eye(N0));
J4 = kron([zeros(Nv-1,1),eye(Nv-1)],eye(N0));

for ii = 1:size(J1,1)
    J11(ii) = find(J1(ii,:)==1);
    J21(ii) = find(J2(ii,:)==1);
end

for ii = 1:size(J3,1)
    J31(ii) = find(J3(ii,:)==1);
    J41(ii) = find(J4(ii,:)==1);
end

% end
for nt = 1:Nt
    for n0 = 1:N0
        J_ntn0((nt-1)*N0+n0,(n0-1)*Nt+nt) = 1;
    end
end
for nr = 1:Nr
    for nt = 1:Nt
        nv = (nr-1)*Nt+nt;
        VP(nv) = TP(nt)+PR(nr);
        %             m = (0:N0-1)*Nt+nt;
        Phi0(nv,:) = exp(1i*2*pi*((nt:Nt:N)-1)*(VP(nv))/2/N/deltaR);
        psi0(nv) = exp(1i*2*pi*(VP(nv))/lambda);
    end
end
for k = 1:K
    t = Tmin + (k-1)*ts + (tc:1/fs:(ts-1/fs));
    for nr = 1:Nr
        Y1(:,k) = diag(exp(-1i*2*pi*(0:N-1)*Tmin/tb))...
            *conj(diag(a(:,k)))*fft(exp(-1i*2*pi*f0*t.').*Yw((nr-1)*N+(1:N),k));
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            Y2((nr-1)*N+(nt-1)*N0+(1:N0),k) = psi0(nv)*diag(Phi0(nv,:))*Y1(nt:Nt:N,k);
            %             Y3((nr-1)*N+(nt:Nt:N),k) = psi0*Phi0*Y1(nt:Nt:N,k);
        end
    end
end
save('data\Preprodata_targetPlane_0409.mat');

YY = reshape(Y2(:,1),[N0,Nt*Nr]);
NFFT = 1024;
for nv = 1:Nt*Nr
    fft_image(:,nv) = fftshift(ifft(YY(:,nv),NFFT));
end
for nf = 1:NFFT
    fft_image2(nf,:) = fftshift(fft(fft_image(nf,:),NFFT));
end
figure
imagesc(abs(fft_image2))
axis xy
Nv = Nr*Nt;
% Z =  unitary_transform_0403(Y2);
Es = Uni_ss_est(Y2,I);

Psi = (Es(J11,:))\(Es(J21,:));
[T1,D] = eig(Psi);
f1 = angle(diag(D))/2/pi/Nt;
ui_est = -f1*N*deltaR;
A = Es*T1;%(:,inds);
for nt = 1:Nt
    for nr = 1:Nr
        nv = (nr-1)*Nt+nt;
        quaterm(nv,:) = (R(1,nr)^2+T(1,nt)^2)/2/OP^2*ui_est.';
        comp0(nv,:) = exp(-1i*2*pi*quaterm(nv,:)/lambda).*exp(-1i*2*pi*(nt-1)*f1.');
    end
end
Comp = kron(comp0,ones(N0,1));

Phi1 = (A(J31,:))\(A(J41,:));
f2 = phase(diag(Phi1))/2/pi;
% [f1,f2] = Uni_ESPRIT_2D(Y2,N0,Nv,I);

vi_est = f2*lambda*OP/dt;

A1 = Comp.*A;
Phi = (A1(J31,:))\(A1(J41,:));
f21 = phase(diag(Phi))/2/pi;
% [f1,f2] = Uni_ESPRIT_2D(Y2,N0,Nv,I);

vi_est1 = f21*lambda*OP/dt;

% ui_est = -angle(phi_est)/2/pi*N*deltaR;
figure(1)
hold on
scatter(ui_est,vi_est,'+')
scatter(ui_est,vi_est1,'*')
% tmp = UniMat(N-1)'*([zeros(N-1,1) eye(N-1)]*permuMat(Nt,N0))*UniMat(N);
% K1 = kron(eye(Nr),real(tmp));
% K2 = kron(eye(Nr),imag(tmp));
% tmp = UniMat(Nv-1)'*([zeros(Nv-1,1) eye(Nv-1)])*UniMat(Nv);
% K3 = kron(real(tmp),eye(N0));
% K4 = kron(imag(tmp),eye(N0));
% omega1 = (K1*Es)\(K2*Es);
% omega2 = (K3*Es)\(K4*Es);
% [V,D] = eig(omega1+1i*omega2);
% ui = -2*atan(real(diag(D)))/2/pi*N*deltaR;
% vi = 2*atan(imag(diag(D)))/2/pi*lambda*R0/dt;
% figure
% scatter(ui,vi)