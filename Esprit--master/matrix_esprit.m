%MATRIX ESPRIT ALOGRITHM
%DOA ESTIMATION BY MATRIX ESPRIT ALOGRITHM
clear all;
close all;
clc;

source_number=2;%��Ԫ��
sensor_number=8;%ԭ��Ԫ��
m=7;%����Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=[pi/4 pi/6].';%�ź�Ƶ��
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%�źŲ���  
d=0.5*l;%��Ԫ���
array_distance=d;%������Ԫ���߷���ļ��
snr=0;%�����

source_doa=[45 60];%�����źŵ�����Ƕ�

A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%��������

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%�����ź�
%��������Z
Z=zeros(m,m);
for i=1:(m-1)
    Z(i+1,i)=1;
end

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�

x1=x(1:m,:);%����1���յ�����
x2=x(2:sensor_number,:);%����2���յ�����

R11=x1*x1'/snapshot_number;%����1������ؾ���
R12=x1*x2'/snapshot_number;%����1��2�Ļ���ؾ���

D1=eig(R11);
u2=mean(D1(1:m-source_number)); %���������ķ���

%ȥ��õ�����1��2���Ժͻ�Э�������
C11=R11-u2*eye(m);
C12=R12-u2*Z;

%��C11��C12���й�������ֵ�ֽ⣬�õ�N��������ֵ
D=eig(C11,C12);
D_big=zeros(1,source_number);
[Y L]=sort(abs(abs(D)-1));
D=D(L);
disp(D);
D_big=[D(2) D(1)];
%����õ��źŵ����
doa =-asin((angle(D_big)/pi))*180/pi %������������






