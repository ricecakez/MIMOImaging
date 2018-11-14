%R_ESPRIT ALOGRITHM
%DOA ESTIMATION BY R_ESPRIT
clear all;
close all;
clc;

source_number=2;%��Ԫ��
sensor_number=8;%��Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=[pi/4 pi/6].';%�ź�Ƶ��
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%�źŲ���  
d=0.5*l;%��Ԫ���
snr=0;%�����

source_doa=[45 60];%�����źŵ�����Ƕ�
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%��������

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�
%�������Qm��Q2L
J=[0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0];%�Ľ��û�����
Qm=(1/sqrt(2))*[eye(4) j*eye(4);J -j*J];
JJ=zeros(snapshot_number,snapshot_number);
for i=1:snapshot_number
      JJ(i,snapshot_number+1-i)=1;
end
Q2L=(1/sqrt(2))*[eye(snapshot_number) j*eye(snapshot_number);JJ -j*JJ];
%����Z����
Jm=[0 0 0 0 0 0 0 1;0 0 0 0 0 0 1 0;0 0 0 0 0 1 0 0;0 0 0 0 1 0 0 0;0 0 0 1 0 0 0 0;0 0 1 0 0 0 0 0;0 1 0 0 0 0 0 0;1 0 0 0 0 0 0 0];
xx=(x').';%������ʸ��ȡ����
Z=[x Jm*xx*JJ];
%����任����Tx
Tx=Qm'*Z*Q2L;
%ͨ���任����Tx�����н�������ת����ʵֵ�ռ�Rt
Rt=Tx*Tx'/(2*snapshot_number);
%��Rt��������ֵ�ֽ⣬�õ����ź��ӿռ�Es
[Ur,Sr,Vr]=svd(Rt);
Es=Ur(:,1:source_number);
%����LS�����ת������
kk=[0 0 0 0 0 0 0].';
K2=[kk eye(sensor_number-1)];%�������K2
%�������Qmm
JJJ=[0 0 1;0 1 0;1 0 0];
Qmm=(1/sqrt(2))*[eye(3) [0 0 0].' j*eye(3);0 0 0 sqrt(2) 0 0 0;JJJ [0 0 0].' -j*JJJ];
%����õ�H1��H2
H1=2*real(Qmm'*K2*Qm);
H2=2*imag(Qmm'*K2*Qm);
%����ת�������M
M=pinv(H1*Es)*(H2*Es);

%�Եõ�����ת���������������ֽ�
[Vm,Dm]=eig(M);
disp(Dm);
beta(1)=-2*atan(Dm(1,1));
beta(2)=-2*atan(Dm(2,2));
disp(beta);
%����õ��źŵ����
doa(1)=real(asin(beta(1)/pi)*180/pi);
doa(2)=real(asin(beta(2)/pi)*180/pi);
disp(doa);
