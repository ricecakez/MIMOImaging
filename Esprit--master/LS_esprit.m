%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALOGRITHM
clear all;
%close all;
clc;

source_number=2;%��Ԫ��
sensor_number=8;%ԭ��Ԫ��
m=7;%����Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=[pi/4 pi/6].';%�ź�Ƶ��
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%�źŲ���  
d=0.5*l;%��Ԫ���

snr=0;%�����/dB
source_doa=[45 60];%�����źŵ�����Ƕ�
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%��������

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�

x1=x(1:m,:);%����1���ܵ�����ʸ��
x2=x(2:(m+1),:);%����2���ܵ�����ʸ��

%�����������ģ�ͽ��кϲ�
X=[x1;x2];
R=X*X'/snapshot_number;

%��R��������ֵ�ֽ�
[U,S,V]=svd(R);
R=R-S(2*m,2*m)*eye(2*m);
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m,:);
Us2=Us((m+1):2*m,:);


%���չ�ʽ�õ���ת�������M
M=pinv(Us1)*Us2;
disp('M');
disp(M);
%�Եõ�����ת���������������ֽ�
[Vm,Dm]=eig(M);
disp(Dm);
Dm=(diag(Dm)).';
doa=-asin(angle(Dm)/pi)*180/pi;
disp(doa);