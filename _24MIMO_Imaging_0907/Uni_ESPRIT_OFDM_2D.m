function [f1,f2,S] = Uni_ESPRIT_OFDM_2D(X,N0,Nt,Nr,K)
%     [Es,~,~] = svds(X,K);
    Es = Uni_ss_est(X,K);
    Nv = Nt*Nr;
    N = N0*Nt;
    J1 = kron(eye(Nv),[eye(N0-1),zeros(N0-1,1)]);
    J2 = kron(eye(Nv),[zeros(N0-1,1),eye(N0-1)]);
    Psi = (J1*Es)\(J2*Es);
    [T,D] = eig(Psi);
    f1 = angle(diag(D))/2/pi/Nt;
%     [f1,inds] = sort(f1);
    A = Es*T;%(:,inds);
    J3 = kron([eye(Nr-1),zeros(Nr-1,1)],eye(N));
    J4 = kron([zeros(Nr-1,1),eye(Nr-1)],eye(N));
    Phi = (J3*A)\(J4*A);
    f2 = angle(diag(Phi))/2/pi/Nt;
    
    A1 = exp(1i*2*pi*(0:N0-1).'*f1.'*Nt);
    A2 = exp(1i*2*pi*(0:Nt-1).'*(f1+f2).');
    A3 = exp(1i*2*pi*(0:Nr-1).'*f2.'*Nt);
    A = krb(A3,krb(A2,A1));
    S = mean(abs(A\X),2);
end