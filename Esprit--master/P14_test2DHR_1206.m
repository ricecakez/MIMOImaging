clear;clc;close all;

K = 20;
L = 20;
F = 3;
omega = [0.48 0.48 0.52]*pi;
nu = [0.48 0.44 0.48]*pi;
C = diag(randn(1,F));
A = exp(1i*(0:K-1).'*omega);
B = exp(1i*(0:L-1).'*nu);
X0 = A*C*B.';
SNR = 30;
% SNR = [0:5:35];
rmse = zeros(size(SNR));
MM = 500;
figure(1)
        scatter(omega,nu,'o')
        hold on
% X = X0;
for nn = 1:length(SNR)
    for mm = 1:MM
        X = awgn(X0,SNR(nn),'measured');
        Y = IEM(K)*conj(X)*IEM(L);
        %% 2D MDF
        % forward-backward averaging
        if mod(K,2) == 0
            K1 = K/2;
            K2 = K+1-K1;
        else
            K1 = (K+1)/2;
            K2 = K1;
        end
        if mod(L,2) == 0
            L1 = L/2;
            L2 = L+1-L1;
        else
            L1 = (L+1)/2;
            L2 = L1;
        end
        X_1 = zeros(K1*L1,K2*L2);
        Y_1 = zeros(K1*L1,K2*L2);
        for k1 = 1:K1
            for k2 = 1:K2
                for l1 = 1:L1
                    for l2 = 1:L2
                        u = (k1-1)*L1+l1;
                        v = (k2-1)*L2+l2;
                        X_1(u,v) = X(k1+k2-1,l1+l2-1);
                        Y_1(u,v) = conj(X(K-(k1+k2-1)+1,L-(l1+l2-1)+1));
                    end
                end
            end
        end
        %         X_ = [];
        %         Y_ = [];
        %         %X_1 = [];
        %         for k2 = 1:K2
        %             X_k2 = [];
        %             Y_k2 = [];
        %             for l1 = 1:L1
        %                 X_l1 = X(:,l1:(L2+l1-1));
        %                 X_k2 = [X_k2;X_l1(k2:K1+k2-1,:)];
        %                 Y_l1 = Y(:,l1:(L2+l1-1));
        %                 Y_k2 = [Y_k2;Y_l1(k2:K1+k2-1,:)];
        %             end
        %             %     X_k22 = kr(B(1:L1,:),A(1:K1,:))*diag(A(k2,:))*C*B(1:L2,:).';
        %             X_ = [X_,X_k2];
        %             Y_ = [Y_,Y_k2];
        %         end
        %         J = [];
        %         Ik1 = eye(K1);
        %         for k1 = 1:K1
        %             J = [J,kron(eye(L1),Ik1(:,k1))];
        %         end
        %         J = J.';
        %         X_1 = J*X_;
        %         Y_1 = J*Y_;
        Z = [X_1;Y_1];
        [U,~,~] = svd(Z);
        U_ = U(:,1:F);
        P = size(U_,1);
        U_1 = U_(1:P/2,:);
        U_2 = U_(P/2+1:end,:);
        [T_1,Phi] = eig(U_1\U_2);
        G = U_1*T_1;
        J3 = kron([eye(K1-1) zeros(K1-1,1)],eye(L1));
        J4 = kron([zeros(K1-1,1) eye(K1-1)],eye(L1));
        
        J1 = kron(eye(K1),[eye(L1-1) zeros(L1-1,1)]);
        J2 = kron(eye(K1),[zeros(L1-1,1) eye(L1-1)]);
        omega_est = angle(diag((J3*G)\(J4*G)));
        nu_est = angle(diag((J1*G)\(J2*G)));
        [omega_est,nu_est] = sort2D(omega_est,nu_est,1);
        rmse(nn) = rmse(nn) + sum((omega_est.' - omega).^2 + (nu_est.' - nu).^2)/2/F;
        figure(1)
        scatter(omega_est,nu_est,'k.')
    end
    rmse(nn) = sqrt(rmse(nn)/MM);
end
figure
semilogy(SNR,rmse,'-o')


