clear;clc;%close all;

c = 3e8;
[Fc,Bw,fs,N,K,df,tb,tc,T,PRF,Nt,Nr,dt,dr] = MIMORadarPara;
lambda = c/Fc;
N_T = fs*T;
N_c = fs*tc;
N0 = N/Nt;      %subcarrier number of interleaved OFDM
[targetRCS,targetPos,targetVel,targetIniPos] = TargetPara;
P = targetIniPos(:);
Qi = targetPos.';
I = size(Qi,2);
T = [(0:Nt-1)*dt;zeros(2,Nt)];
R = [(0:Nr-1)*dr;zeros(2,Nr)];
Nv = Nt*Nr;
R0 = rangeangle(P);
u0 = P/R0;
qi = Qi-P;
for nt = 1:Nt
    for nr = 1:Nr
        nv = nt + (nr-1)*Nt;
        Rc(nv) = -2*R0 + P(1)/R0*dt*(nv-1)...
            -1/2*(P(2)^2+P(3)^2)/R0^3*dt^2*((nt-1)^2+((nr-1)*Nt)^2);%...
         Rc1(nv) = (rangeangle(T(:,nt),P)+rangeangle(R(:,nr),P));
%             +(P(2)^2+P(3)^2)/R0^3*dt^2*(nv-1)^2;
%         R_vc(nv) = -P(1)/R0*dt*(nv-1)+1/2*(1/R0-P(1)^2/R0^3)*dt^2*(nv-1)^2;
%         R_vc1(nv) = -P(1)/R0*dt*(nv-1);
        
    end
end
figure
plot(Rc - Rc1)
d0 = c/2/Bw;
Rn = (-N/2:(N/2-1))/N*d0*N0;
dR1 = Rn(2)-Rn(1);
% tmp4 = R_vc - R_vc1;
% figure(2)
% plot(R_vc)
% hold on
% plot(R_vc1)
for i = 1:I
    Ri = rangeangle(Qi(:,i));
    at = Qi(1,i)*dt/Ri;
    bt = 1/2*(1/Ri-Qi(1,i)^2/Ri^3)*dt^2;
    ar = Qi(1,i)*dr/Ri;
    br = 1/2*(1/Ri-Qi(1,i)^2/Ri^3)*dr^2;
    TQi_est = Ri - at*(0:Nt-1) + bt*(0:Nt-1).^2;
    TQi = rangeangle(T,Qi(:,i));
    tmp = max(abs(TQi_est - TQi));
    QiR_est = Ri - ar*(0:Nr-1) + br*(0:Nr-1).^2;
    QiR = rangeangle(R,Qi(:,i));
    tmp1 = max(abs(QiR_est - QiR));
    Qi1 = [P(1) + (Qi(1,i)-P(1))/2;Qi(2,i);Qi(3,i)];
    Ri1 = rangeangle(Qi1);
    tmp2 = Ri1 - Ri;
    vi(i) = (qi(:,i) - (Ri-R0)*u0).'*[1;0;0];
    vi1(i) = (qi(:,i) - (qi(:,i).'*u0)*u0).'*[1;0;0];
    %     tmpe = abs((1/R0-P(1)^2/R0^3)-(1/Ri-Qi(1,i)^2/Ri^3))*dt^2
    for nr = 1:Nr
        for nt = 1:Nt
            nv = (nr-1)*Nt+nt;
            V(:,nv) =  [dt*(nv-1);0;0];
            
            %             TQiR_est1(i,nv) = 2*Ri - at*(nv-1) + bt*(nv-1)^2;
            %             TQiR_est(i,nv) = TQi_est(nt) + QiR_est(nr);
            TQiR(i,nv) = TQi(nt) + QiR(nr);
            R0v(nv) = rangeangle(V(:,nv),P);
            R1(i,nv) = -TQiR(i,nv)+Rc1(nv);
%             nn(i,nv) = round((-(R1(i,nv)/2)-Rn(1))/dR1)+1;%v(nv);
%             R1_est(i,nv) = 2*(Ri-R0-((Qi(1,i)-P(1)))/2/Ri1*dt*(nv-1));
%             nn_est(i,nv) = round((-R1_est(i,nv)/2-Rn(1))/dR1)+1;
%             R1_est(i,nv) = 2*((Ri-R0)-((Qi(1,i)-P(1))/2)*dt/Ri*(nv-1));%+1/2*(Qi(2,i)^2+Qi(3,i)^2)/Ri1^3*dt^2*(nv-1)^2);
%             R1_est2(i,nv) = 2*((Ri-R0)-((Qi(1,i)-P(1))/2)*dt/Ri*(nv-1));%+1/2*(P(2)^2+P(3)^2)/R0^3*dt^2*(nv-1)^2); 
%            
            R1_est(i,nv) = 2*(R0-Ri) + (nv-1)*dt/R0*vi1(i)...
                +(qi(:,i).'*u0)*dt^2/2/R0^2*((nt-1)^2+((nr-1)*Nt)^2);
             tmp3(i,nv) = R1(i,nv) - R1_est(i,nv);
            R1_est1(i,nv) = 2*(R0-Ri) + (nv-1)*dt/R0*vi(i)...%;
                + (Ri-R0)*dt^2/2/R0^2*((nt-1)^2+((nr-1)*Nt)^2);
            tmp4(i,nv) = R1(i,nv) - R1_est1(i,nv);
%             tmpe(i,nv) = P(1)*(Ri/R0-1)/Ri*dt*(nv-1);
%             delta_Ri(i) = Ri-R0;
            %             tmp(i,nv) = abs(TQiR(i,nv) - TQiR_est(i,nv));
            %             tmp1(i,nv) = abs(TQiR(i,nv) - TQiR_est1(i,nv));
        end
    end
    figure(3)
    plot((0:Nv-1),R1(i,:),(0:Nv-1),R1_est(i,:),(0:Nv-1),R1_est1(i,:))
    hold on
    figure(4)
    plot(tmp3(i,:))
    hold on
    plot(tmp4(i,:))
end
% % x = sum(exp(-1i*2*pi*(R1_est-2*R_vc-2*R0)/lambda),1);
% % x1 = [x,zeros(1,Nv*3)];
% % plot(abs(fftshift(ifft(x,Nv*4))));
% % fx = fftshift(frft(x1,-0.01));
% % figure
% % plot(abs(fx))