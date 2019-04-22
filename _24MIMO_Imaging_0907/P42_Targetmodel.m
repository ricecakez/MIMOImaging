clear;clc;close all;

c = 3e8;
SNR = 5;
[f0,lambda,B,fs,N,K,df,tb,tc,ts,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
deltaR = c/2/B;
[rho,SP,VT,P] = TargetPara;
% SP = SP(1,:);
I = size(SP,1);
Q = SP.';
P = P.';
q = Q - P;
SP1 = P + [q(1,:)/2;q(2,:);q(3,:)];
Nr = 16;
T = [dt;0;0]*(0:Nt-1);
R = [dr;0;0]*(0:Nr-1);
[OP,theta0] = rangeangle(P);
OQ = rangeangle(Q,T(:,1));
Ri1 = rangeangle(SP1,T(:,1));
u0 = P/OP;
u01 = [cosd(theta0(1))*cosd(theta0(2));sind(theta0(1))*cosd(theta0(2));...
    sind(theta0(2))];
for i = 1:I
    close all;
% i = 2;
ui(i) = q(:,i).'*u0;
vi(i) = (q(:,i)-u0*ui(i)).'*[1;0;0];
% (vi*dt/R0)/N/deltaR
% ui/N/deltaR
TP = rangeangle(T,P);
PR = rangeangle(R,P);
TQ = rangeangle(T,Q(:,i));
QR = rangeangle(R,Q(:,i));
OQ = rangeangle(Q(:,i));
R1 = (TP - OP) + (OQ - TQ);
R2 = (PR - OP) + (OQ - QR);
RR = kron(ones(1,Nr),R1) + kron(R2,ones(1,Nt));
for nt = 1:Nt
%     R2(nt) = 2*P.'*T(:,nt)*((TP(nt)+OP)-(TQ(nt)+OQ))...
%         /((TP(nt)+OP)*(TQ(nt)+OQ))...
%         +T(:,nt).'*T(:,nt)*(-(TP(nt)+OP)+(TQ(nt)+OQ))...
%         /((TP(nt)+OP)*(TQ(nt)+OQ))...
%         +2*Q(:,1).'*T(:,nt)/(TQ(nt)+OQ);
%     R3(nt) = -u0.'*T(:,nt)*Q(:,1).'*u0/OP+T(:,nt).'*T(:,nt)*Q(:,1).'*u0/2/OP^2 ...
%         +Q(:,1).'*T(:,nt)/OP;
    R3(nt) = vi(i)*dt*(nt-1)/OP+T(:,nt).'*T(:,nt)*ui(i)/2/OP^2;
    
end
for nr = 1:Nr
    R4(nr) = (q(:,i)-q(:,i).'*u0*u0).'*R(:,nr)/OP+R(:,nr).'*R(:,nr)*ui(i)/2/OP^2 ;
end
RR1 = kron(ones(1,Nr),R3) + kron(R4,ones(1,Nt));
figure
ss1 = R1(2:end)- R1(1:end-1);

plot(ss1(2:end) - ss1(1:end-1))
ss3 = R3(2:end)-R3(1:end-1);
hold on
plot(ss3(2:end) - ss3(1:end-1))
% figure
ss2 = (R2(2:end)-R2(1:end-1));
figure

plot(ss2(2:end) - ss2(1:end-1))
% hold on
ss4 = (R4(2:end)-R4(1:end-1));
hold on
plot(ss4(2:end) - ss4(1:end-1))
% figure
rr(i) = mean(RR(2:end)-RR(1:end-1))
% hold on
rr1(i) = mean(RR1(2:end)-RR1(1:end-1))

% figure
% plot(R2(2:end)-R2(1:end-1))
% hold on
% plot(R4(2:end)-R4(1:end-1))
% max(abs(R2-R1))
nn(i) = max(abs(RR1-RR));
end
% max(abs(R4-R2))
R_PR = rangeangle(R,P);
% plot([R_TP,R_PR])
for nr = 1:Nr
    for nt = 1:Nt
        nv = (nr-1)*Nt+nt;
%         V(:,nv) = [dt;0;0]*(nv-1);
        R_i(:,nv) = (TP(nt)+R_PR(nr))-(rangeangle(Q,T(:,nt))+rangeangle(Q,R(:,nr)));
        R_i1(:,nv) = vi*(T(1,nt)+R(1,nr))/OP - 2*ui.'+ ui.'*(T(1,nt)^2+R(1,nr)^2)/OP^2/2;
%         R_TQR1(:,nv) = 2*Ri1-SP1(1,:)./Ri1*(T(1,nt)+R(1,nr))...
%             +1/2*(ones(1,I)./Ri1-(SP1(1,:).^2)./(Ri1.^3))*(T(1,nt)^2+R(1,nr)^2);
%         R_adjust(:,nv) = -P(1)/R0*(V(1,nv))...
%             +1/2*(1/R0-P(1)^2/R0^3)*(T(1,nt)^2+R(1,nr)^2);%...
%             %-(1/R0-P(1)^2/R0^3)*(T(1,nt)+R(1,nr))^2;
%         R_i(:,nv) = R_TQR1(:,nv) - R_adjust(:,nv);
%         R_VQ2(:,nv) = 2*(Ri1-(SP1(1,:)./Ri1)*V(1,nv)...
%             +1/2*(ones(1,I)./Ri1-(SP1(1,:).^2)./(Ri1.^3))*(V(1,nv))^2);
%         R_VQ(:,nv) = 2*rangeangle(SP1,V(:,nv));
        
%         R_modi(nv) = 2*(-P(1)/R0*V(1,nv));%+1/2*(1/R0-P(1)^2/R0^3)*V(1,nv)^2);
%         R_i(:,nv) = R_VQ1(:,nv) - R_modi(nv);
%         R_i1(:,nv) = 2*(Ri1-Q(1,:)./Ri1/2*V(1,nv));
%         R_TPR(nv) = R_TP(nt) + R_PR(nr);
%         R_VP(nv) = 2*rangeangle(V(:,nv),P);
%         eps1(:,nv) = R_TPR(nv) - R_TQR(:,nv);
%         eps2(:,nv) = R_VP(nv) - R_VQ(:,nv);
%         eps3(:,nv) = vi*(T(1,nt)+R(1,nr))/R0-2*ui+ui*(T(1,nt)^2+R(1,nr)^2)/2/R0^2;
%         eps4(:,nv) = 2*(vi*V(1,nv)/R0-ui+ui*(V(1,nv)^2)/2/R0^2);
    end
end
% figure
% plot(R_TQR(1,2:end)-R_TQR(1,1:end-1))
% hold on
% plot(R_TQR1(1,2:end)-R_TQR1(1,1:end-1))

% figure
% plot(R_VQ1(1,2:end)-R_VQ1(1,1:end-1))
% hold on
% % plot(R_VQ2(1,2:end)-R_VQ2(1,1:end-1))
% plot(R_VQ(1,2:end)-R_VQ(1,1:end-1))
for i = 1:I
figure
plot((R_i(i,2:end)-R_i(i,1:end-1)))
hold on
plot((R_i1(i,2:end)-R_i1(i,1:end-1)),'--')
end
% figure
% plot(R_i - R_i1)