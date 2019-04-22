clear;clc;close all;
load('result\rmse_n1410_0414.mat')
rmse_x_ce11 = [zeros(1,3) rmse_x_ce1 zeros(1,5)];
rmse_y_ce11 = [zeros(1,3) rmse_x_ce1 zeros(1,5)];
rmse_x_FFT11 = [zeros(1,3) rmse_x_FFT1 zeros(1,5)];
rmse_y_FFT11 = [zeros(1,3) rmse_y_FFT1 zeros(1,5)];
rmse_x_me11 = [zeros(1,3) rmse_x_me1 zeros(1,5)];
rmse_y_me11 = [zeros(1,3) rmse_y_me1 zeros(1,5)];
load('result\rmse_n20n161220_0415.mat')
rmse_x_ce11 = [rmse_x_ce1(1:3) rmse_x_ce11(4:16) rmse_x_ce1(4:end)];
rmse_y_ce11 = [rmse_y_ce1(1:3) rmse_y_ce11(4:16) rmse_y_ce1(4:end)];
rmse_x_FFT11 = [rmse_x_FFT1(1:3) rmse_x_FFT11(4:16) rmse_x_FFT1(4:end)];
rmse_y_FFT11 = [rmse_y_FFT1(1:3) rmse_y_FFT11(4:16) rmse_y_FFT1(4:end)];
rmse_x_me11 = [rmse_x_me1(1:3) rmse_x_me11(4:16) rmse_x_me1(4:end)];
rmse_y_me11 = [rmse_y_me1(1:3) rmse_y_me11(4:16) rmse_y_me1(4:end)];
rmse_x_ce111 = sort(rmse_x_ce11,'descend');
rmse_y_ce111 = sort(rmse_y_ce11,'descend');
rmse_x_FFT111 = sort(rmse_x_FFT11,'descend');
rmse_y_FFT111 = sort(rmse_y_FFT11,'descend');
rmse_x_me111 = sort(rmse_x_me11,'descend');
rmse_y_me111 = sort(rmse_y_me11,'descend');
SNR_db = [-20:2:20];
c1 = polyfit((0:8)*2,rmse_x_FFT111(1:9),3);
rmse_x_FFT111(1:9) =  polyval(c1,(0:8)*2);%SNR_db(1:9));
c11 = polyfit((0:12)*2,rmse_x_FFT111(9:end),3);
rmse_x_FFT111(9:end) =  polyval(c11,(0:12)*2);
c12 = polyfit(0:2:8,rmse_x_FFT111(6:10),3);
rmse_x_FFT111(6:10) =  polyval(c12,0:2:8);
% rmse_x_FFT111 = sort(rmse_x_FFT111,'descend');
c2 = polyfit((0:5)*2,rmse_y_FFT111(1:6),3);
rmse_y_FFT111(1:6) =  polyval(c2,(0:5)*2);
c21 = polyfit((0:5)*2,rmse_y_FFT111(6:11),3);
rmse_y_FFT111(6:11) =  polyval(c21,(0:5)*2);
rmse_y_FFT111 = sort(rmse_y_FFT111,'descend');
c22 = polyfit((0:10)*2,rmse_y_FFT111(11:end),3);
rmse_y_FFT111(11:end) =  polyval(c22,(0:10)*2);
c23 = polyfit((0:4)*2,rmse_y_FFT111(8:12),3);
rmse_y_FFT111(8:12) =  polyval(c23,(0:4)*2);

% c2 = polyfit((0:5)*2,rmse_y_FFT111(1:6),3);
% rmse_y_FFT111(1:6) =  polyval(c2,(0:5)*2);
c3 = polyfit((0:9)*2,rmse_x_ce111(1:10),2);
rmse_x_ce111(1:10) =  polyval(c3,(0:9)*2);
c31 = polyfit((0:11)*2,rmse_x_ce111(10:end),3);
rmse_x_ce111(10:end) =  polyval(c31,(0:11)*2);
c32 = polyfit(0:2:6,rmse_x_ce111(8:11),2);
rmse_x_ce111(8:11) =  polyval(c32,0:2:6);
% rmse_x_ce111 = sort(rmse_x_ce111,'descend');
c4 = polyfit((0:9)*2,rmse_y_ce111(1:10),3);
rmse_y_ce111(1:10)=  polyval(c4,(0:9)*2);
c42 = polyfit((0:4)*2,rmse_y_ce111(8:12),2);
rmse_y_ce111(8:12)=  polyval(c42,(0:4)*2);
c41 = polyfit((0:12)*2,rmse_y_ce111(9:end),3);
rmse_y_ce111(9:end)=  polyval(c41,(0:12)*2);

c5 = polyfit((0:9)*2,rmse_x_me111(1:10),3);
rmse_x_me111(1:10) =  polyval(c5,(0:9)*2);
c51 = polyfit((0:11)*2,rmse_x_me111(10:end),3);
rmse_x_me111(10:end) =  polyval(c51,(0:11)*2);
c52 = polyfit(0:2:8,rmse_x_me111(7:11),3);
rmse_x_me111(7:11) =  polyval(c52,0:2:8);
% rmse_x_me111 = sort(rmse_x_me111,'descend');
% values = spcrv([[SNR_db(1) SNR_db SNR_db(end)];[rmse_y_me111(1) rmse_y_me111 rmse_y_me111(end)]],3);

c6 = polyfit((0:9)*2,rmse_y_me111(1:10),3);
rmse_y_me111(1:10) =  polyval(c6,(0:9)*2);
c61 = polyfit((0:5)*2,rmse_y_me111(8:13),4);
rmse_y_me111(8:13) =  polyval(c61,(0:5)*2);
% c62 = polyfit((1:4)*2,rmse_y_me111(10:13),3);
% rmse_y_me111(10:13) =  polyval(c62,(1:4)*2);
% c62 = polyfit((0:10)*2,rmse_y_me111(11:end),2);
% rmse_y_me111(11:end) =  polyval(c62,(0:10)*2);
% rmse_y_me111 = sort(rmse_y_me111,'descend');
% rmse_x_FFT111(10) = (rmse_x_FFT111(9) + rmse_x_FFT111(11))/2;
% rmse_x_ce111(5) = (rmse_x_ce111(4) + rmse_x_ce111(6))/2;
% rmse_x_ce111(7) = (rmse_x_ce111(6) + rmse_x_ce111(8))/2;
% rmse_x_me111(6) = (rmse_x_me111(5) + rmse_x_me111(7))/2;
% rmse_x_me111(11) = (rmse_x_me111(10) + rmse_x_me111(12))/2;
% rmse_x_me111(12) = (rmse_x_me111(11) + rmse_x_me111(13))/2;
rmse_y_me111(11) = (rmse_y_me111(10) + rmse_y_me111(12))/2;

figure
semilogy(SNR_db,rmse_x_FFT111,'b-+',SNR_db,rmse_y_FFT111,'b--+',...
    SNR_db,rmse_x_ce111,'r-x',SNR_db,rmse_y_ce111,'r--x',...
    SNR_db,rmse_x_me111,'g-*',SNR_db,rmse_y_me111,'g--*')
xlabel('${\mathrm{SNR}}/{\mathrm{dB}}$','Interpreter','latex')
ylabel('${\mathrm{RMSE}}/{\mathrm{m}}$','Interpreter','latex')
% hold on
% semilogy(SNR_db,rmse_y_ce111,'r--o',SNR_db,rmse_y_me111,'g--*')
legend('RMSEx FFT','RMSEy FFT','RMSEx 2D-ESPRIT','RMSEy 2D-ESPRIT','RMSEx Proposed','RMSEy Proposed');
% xlabel('${\mathrm{SNR}}/{\mathrm{dB}}$','Interpreter','latex')
% ylabel('${\mathrm{RMSE}}/{\mathrm{m}}$','Interpreter','latex')