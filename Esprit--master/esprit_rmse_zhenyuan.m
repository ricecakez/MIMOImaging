%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT
clear all;
%close all;
clc;

bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%��Ԫ��
m=[2 3 5 7 9 11 13 15 17 19];%����Ԫ��
aaa=zeros(1,300);
for k=1:300

source_number=1;%��Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=pi/4;%�ź�Ƶ��
l=2*pi*3e8/w;%�źŲ���  
d=0.5*l;%��Ԫ���
snr=-2;%�����(dB)
source_doa=50;%�źŵ�����Ƕ�
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';

s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%���˸�˹������������н����ź�

x1=x(1:m(kk),:);%����1���ܵ�����ʸ��
x2=x(2:(m(kk)+1),:);%����2���ܵ�����ʸ��

%�����������ģ�ͽ��кϲ�
X=[x1;x2];
R=X*X'/snapshot_number;
%��R��������ֵ�ֽ�
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m(kk),:);
Us2=Us((m(kk)+1):2*m(kk),:);

%���չ�ʽ�õ���ת�������M
M=pinv(Us1)*Us2;
disp('M');
disp(M);
%�Եõ�����ת���������������ֽ�
[Vm,Dm]=eig(M);
disp(Dm);
estimated_source_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_source_doa);

aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%�����������ͱ�׼ƫ��
%E_source_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
E_source_doa=mean(aaa);
disp(E_source_doa);

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

plot(sensor_number,bbb(1,:),'k*-');
save LS_ESPRIT_zhenyuan_RMSE.mat;
%TLS_ESPRIT
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%��Ԫ��
m=[2 3 5 7 9 11 13 15 17 19];%����Ԫ��

aaa=zeros(1,300);
for k=1:300
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';

s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%���˸�˹������������н����ź�
x1=x(1:m(kk),:);%����1���ܵ�����ʸ��
x2=x(2:(m(kk)+1),:);%����2���ܵ�����ʸ��

%�����������ģ�ͽ��кϲ�
X=[x1;x2];
R=X*X'/snapshot_number;
%��R��������ֵ�ֽ�
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m(kk),:);
Us2=Us((m(kk)+1):2*m(kk),:);
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
estimated_source_doa=-asin(angle(Dm(1,1))/pi)*180/pi;

aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%�����������ͱ�׼ƫ��
E_source_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
disp('E_source_doa');

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);
hold on
plot(sensor_number,bbb(1,:),'rs-');
save TLS_ESPRIT_zhenyuan_rmse.mat;
%TAM ALOGRITHM
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%��Ԫ��

aaa=zeros(1,300);
for k=1:300
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%���˸�˹������������н����ź�

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
B1=B(1:(sensor_number(kk)-1),:);
B2=B(2:sensor_number(kk),:);
%��ȡUs������������Us1��Us2
%Us1=Us(1:(sensor_number(kk)-1),:);
%Us2=Us(2:sensor_number(kk),:);
%�������Ϲ�ϵ�õ���С���˽�
%D=pinv(Us1*((Ss)^(1/2)))*Us2*((Ss)^(1/2));
D=pinv(B1)*B2;
%��D���������ֽ⣬������ֵ�ɵõ���Ӧ��N���źŵĵ����
[Vd,Dd]=eig(D);
disp(Dd);
estimated_source_doa=-asin(angle(Dd(1,1))/pi)*180/pi;
%�����źŵ��﷽���
aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%�����������ͱ�׼ƫ��
E_source_doa=sum(aaa(1,:))/300;%��300������ľ�ֵ
disp('E_source_doa');

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%��300������ľ��������
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

hold on
plot(sensor_number,bbb(1,:),'bd-');

save TAM_zhenyuan_RMSE.mat;

legend('LS-ESPRIT','TLS-ESPRIT','TAM');
xlabel('��Ԫ��Ŀ');
ylabel('���ƾ��������');
grid on;





