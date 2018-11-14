clear all;clc;close all;
format long %��������ʾΪ�����Ϳ�ѧ����
N=300;%������
doa=[10 30 60]/180*pi; %�źŵ����
w=[pi/6 pi/4 pi/8]';%�ź�Ƶ��
M=8;%��Ԫ��
P=length(w); %�ź�Դ����
c=3e8;       %����
lambda=c/w;  %����
d=lambda/2;  %��Ԫ���
snr=0;
B=zeros(P,M); %����һ��P��M�е�0����
for i=1:P
B(i,:)=exp(-j*2*pi*d*sin(doa(i))/lambda*[0:M-1]); %����ֵ A
end
B=B.';
S=2*exp(j*(w*[1:N])); %�����ź�
x=B*S;  
x=x+awgn(x,snr);%�����˹������
R =x*x'/N;%����غ���
Rr=real(R);%ʵ��
Ri=imag(R);%�鲿
N=M;
J=zeros(N,N);
for k=1:N
    J(k,N+1-k)=1;
end
J=J;
T=Rr+J*Ri;  %�����µ�ʵ�ϳɾ���
[V,D]=eig(T);  %����ֵ�ֽ�
Us=V(:,M-P+1:M);
U1=Us(1:M-1,:);
U2=Us(2:M,:);% ���������ųɵ������ӿռ�
[p,q]=eig(pinv(U1)*U2); % ������С���˷� �����ת�����ϵ����Ȼ����������ֽ�
for i=1:P;
    alpha(i)=real(asin(j*(log(q(i,i)))*lambda/(-2*pi*d))*180/pi);
end
alpha
stem(alpha,ones(1,P),'filled');grid;
axis([-90 90 0 2]);
text(alpha(3)-4,1.1,num2str(alpha(3)));text(alpha(3)-5,1.3,'�ź�1');
text(alpha(2)-4,1.1,num2str(alpha(2)));text(alpha(2)-5,1.3,'�ź�2');
text(alpha(1)-4,1.1,num2str(alpha(1)));text(alpha(1)-5,1.3,'�ź�3');
ylabel('DOA���ƵĽǶ�ֵ  ');
xlabel('�Ƕ�');
title('ESPRIT�㷨DOA����');