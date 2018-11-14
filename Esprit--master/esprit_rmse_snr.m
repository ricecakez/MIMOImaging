%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALOGRITHM
clear all;
%close all;
clc;

bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%�����(dB)

aaa=zeros(1,300);
for k=1:300

source_number=1;%��Ԫ��
sensor_number=8;%ԭ��Ԫ��
m=7;%����Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=pi/4;%�ź�Ƶ��
l=2*pi*3e8/w;%�źŲ��� 
d=0.5*l;%��Ԫ���

source_doa=50;%�źŵ�����Ƕ�

A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^((snr(kk)/2)/10)*exp(j*w*[0:N_x-1]);%�����ź�

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�

x1=x(1:m,:);%����1���ܵ�����ʸ��
x2=x(2:(m+1),:);%����2���ܵ�����ʸ��

%�����������ģ�ͽ��кϲ�
X=[x1;x2];
R=X*X'/snapshot_number;
%��R��������ֵ�ֽ�
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
estimated_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
plot(snr,bbb(1,:),'k*-');
save LS_ESPRIT_snr_rmse.mat;
%TLS_ESPRIT
bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%�����(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kk)/10))*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�

x1=x(1:m,:);%����1���ܵ�����ʸ��
x2=x(2:(m+1),:);%����2���ܵ�����ʸ��

%�����������ģ�ͽ��кϲ�
X=[x1;x2];
R=X*X'/snapshot_number;
%��R��������ֵ�ֽ�
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m,:);
Us2=Us((m+1):2*m,:);
Us12=[Us1 Us2];
%�γɾ���Us12
Us12=[Us1,Us2];
%�ԡ�Us12'*Us12�����������ֽ⣬�õ�����E
[E,Sa,Va]=svd(Us12'*Us12);
disp('E');
disp(E);
disp(Sa);
%��E�ֽ�Ϊ�ĸ�С����
E11=E(1,1);
E12=E(1,2);
E21=E(2,1);
E22=E(2,2);
%���չ�ʽ�õ���ת�������M
M=-(E12*(inv(E22)));
disp('M');
disp(M);
%�Եõ�����ת���������������ֽ�
[Vm,Dm]=eig(M);
disp(Dm);
estimated_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'rs-');
save TLS_ESPRIT_snr_rmse.mat;
%TAM
bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%�����(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kk)/10))*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�

R=x*x'/snapshot_number;
disp(R);
%�Խ��յ�������Э��������������ֵ�ֽ�,�ֽ�õ��ź��ӿռ估������ֵ
[U,S,V]=svd(R);
Us=U(:,1:source_number);
Ss=S(1:source_number,1:source_number);
disp(Ss);
Vs=V(:,1:source_number);
%������ת�����ӿռ�˼�빹�����B
B=Us*(Ss^(1/2));
B1=B(1:(sensor_number-1),:);
B2=B(2:sensor_number,:);
%��ȡUs������������Us1��Us2
%Us1=Us(1:(sensor_number-1),:);
%Us2=Us(2:sensor_number,:);
%�������Ϲ�ϵ�õ���С���˽�
%D=pinv(Us1*((Ss)^(1/2)))*Us2*((Ss)^(1/2));
D=pinv(B1)*B2;
%��D���������ֽ⣬������ֵ�ɵõ���Ӧ��N���źŵĵ����
[Vd,Dd]=eig(D);
disp(Dd);
estimated_doa=-asin(angle(Dd(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'bd-');
save TAM_snr_rmse.mat
%R_ESPRIT
bbb=zeros(1,11);
for kkk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%�����(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kkk)/10))*exp(j*w*[0:N_x-1]);%�����ź�

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�
%�������Qm��Q2L
J=[0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0];%�Ľ��û�����
Qm=(1/sqrt(2))*[eye(4) j*eye(4);J -j*J];
JJ=zeros(snapshot_number,snapshot_number);
for ii=1:snapshot_number
      JJ(ii,snapshot_number+1-ii)=1;
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
%����õ��źŵ����
estimated_doa =real(asin(beta(1)/pi)*180/pi); %������������
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kkk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'gx-');
save resprit_snr_rmse.mat
%RB_ESPRIT
bbb=zeros(1,11);
for kkk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%�����(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kkk)/10))*exp(j*w*[0:N_x-1]);%�����ź�

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�
%����Ȩ����W
%ac=zeros(sensor_number,sensor_number);
%for k=0:(sensor_number-1)
%gama(k+1)=k*2*pi/sensor_number;
%ac(:,k+1)=exp(-j*((sensor_number-1)/2)*(gama(k+1)/pi))*exp(j*(0:sensor_number-1)'*(gama(k+1)/pi));
%end
%W=ac;
W=[exp(-j*(sensor_number-1)/2*2*0*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*0*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*1*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*1*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*2*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*2*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*3*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*3*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*4*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*4*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*5*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*5*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*6*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*6*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*7*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*7*pi/sensor_number)];
disp(W);
%�����н��յ����ݴӸ���ת����ʵ��
Y=W'*x;
Y1=[real(Y) imag(Y)];
disp(Y1);
R1=Y1*Y1'/(2*snapshot_number);
[U,S,V]=svd(R1);
disp(U);
Es=U(:,1:source_number);
disp(Es);
KS1=[1 cos(pi/sensor_number) cos(2*pi/sensor_number) cos(3*pi/sensor_number) cos(4*pi/sensor_number) cos(5*pi/sensor_number) cos(6*pi/sensor_number) cos(7*pi/sensor_number)];
KS1=diag(KS1);
for i=1:sensor_number-1
    KS1(i,i+1)=cos(i*pi/sensor_number);
end
KS2=[0 sin(pi/sensor_number) sin(2*pi/sensor_number) sin(3*pi/sensor_number) sin(4*pi/sensor_number) sin(5*pi/sensor_number) sin(6*pi/sensor_number) sin(7*pi/sensor_number)];
KS2=diag(KS2);
for ii=1:sensor_number-1
    KS2(ii,ii+1)=sin(ii*pi/sensor_number);
end

M=pinv(KS1*Es)*(KS2*Es);
%�Եõ�����ת���������������ֽ�
[Vm,Dm]=eig(M);
disp(Dm);
beta(1)=-2*atan(Dm(1,1));
%����õ��źŵ����
estimated_doa =real(asin(beta(1)/pi)*180/pi); %������������
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kkk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'mo-');
save RB_esprit_snr_rmse.mat

legend('LS-ESPRIT','TLS-ESPRIT','TAM','ʵֵ�ռ�ESPRIT','ʵֵ�����ռ�ESPRIT');
xlabel('����ȣ�snr��/dB');
ylabel('���ƾ��������');
grid on;


