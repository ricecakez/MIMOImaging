clear variables;
close all;  clc;

%�����(dB)
snr=-10:2:10;
lsnr = length(snr);
RMSE=zeros(1,lsnr);
test_num = 300;

source_number=1; % ��Դ��
sensor_number=8; % ԭ��Ԫ��
m=7;             % ����Ԫ��
N_x=1024;        % �źų���
snapshot_number=N_x; % ������
w=pi/4;              % �ź�Ƶ��
l=2*pi*3e8/w;        % �źŲ���  
d=0.5*l;             % ��Ԫ���
array_distance=d;    % ������Ԫ���߷���ļ��

source_a= 50;  % �����źŵ�����Ƕ�
source_doa= source_a *pi/180;  % �����źŵ�����Ƕ�
% ��������
A=[exp(-1i*(0:sensor_number-1)*d*2*pi*sin(source_doa)/l)].';

for isnr = 1:lsnr
    RMSE_i = zeros(1, test_num);
    for iter = 1:test_num
        % �����ź�
        s=sqrt(10.^(snr(isnr)/10))*exp(1i*w*[0:N_x-1]); 
        %x=awgn(s,snr);
        % ���˸�˹������������н����ź�
        x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+1i*randn(sensor_number,N_x));
        x1=x(1:m,:);     % ����1���ܵ�����ʸ��
        x2=x(2:(m+1),:); % ����2���ܵ�����ʸ��

        % �����������ģ�ͽ��кϲ�
        Rxx = x1*x1'/snapshot_number;
        Rxy = x1*x2'/snapshot_number;
        [V,D] = svd(Rxx);
        ev = diag(D);emin=ev(end);
        Cxx = Rxx - emin*eye(m);
        z=[zeros(m-1,1),eye(m-1);0,zeros(1,m-1)];
        Cxy=Rxy-emin*z;
        [V,D] = eig(Cxx,Cxy);
        z = diag(D);
        
        
%         X = [x1; x2];
%         R=X*X'/snapshot_number;
%         % ��R��������ֵ�ֽ�
%         [U,S,V]=svd(R);
%         R=R-S(2*m,2*m)*eye(2*m);
%         [U,S,V]=svd(R);
%         Us=U(:,1:source_number);
%         % disp(Us);
%         Us1=Us(1:m,:);
%         Us2=Us((m+1):2*m,:);
%         %Us12=[Us1 Us2];
%         % �γɾ���Us12
%         Us12=[Us1,Us2];
%         % �ԡ�Us12'*Us12�����������ֽ⣬�õ�����E
%         [E,Sa,Va]=svd(Us12'*Us12);
%         % disp('E');
%         % disp(E);
%         % disp(Sa);
%         % �� E_4x4 �ֽ�Ϊ�ĸ�С����
%         E11=E(1,1);
%         E12=E(1,2);
%         E21=E(2,1);
%         E22=E(2,2);
%         % ���չ�ʽ�õ���ת�������M
%         M=-(E12*(inv(E22)));
%         % disp('M');
%         % disp(M);
%         % �Եõ�����ת���������������ֽ�
%         [Vm,Dm]=eig(M);
%         % disp(Dm);
%         Dm=(diag(Dm)).';
        doa=-asin(angle(z(1))/pi)*180/pi;
        disp(doa);
        RMSE_i(iter)=sqrt((doa-source_a)^2);
    end
    RMSE(isnr) = sum(RMSE_i)/test_num;
end

%% plot 
figure;
plot(snr,RMSE,'m-s');
%legend('elevation','azimuth');
xlabel('SNR (dB)');
ylabel('RMSE of(degree)');
%axis([-100 100 -40 60]);
%set(gca, 'XTick',[-100:20:100]);
grid on;  