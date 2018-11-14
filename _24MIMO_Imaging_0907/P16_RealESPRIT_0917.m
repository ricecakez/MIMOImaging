clear all;clc;close all;
format long %将数据显示为长整型科学计数
N=300;%快拍数
doa=[10 30 60]/180*pi; %信号到达角
w=[pi/6 pi/4 pi/8]';%信号频率
M=8;%阵元数
P=length(w); %信号源个数
c=3e8;       %光速
lambda=c/w;  %波长
d=lambda/2;  %阵元间距
snr=0;
B=zeros(P,M); %创建一个P行M列的0矩阵
for i=1:P
B(i,:)=exp(-j*2*pi*d*sin(doa(i))/lambda*[0:M-1]); %矩阵赋值 A
end
B=B.';
S=2*exp(j*(w*[1:N])); %仿真信号
x=B*S;  
x=x+awgn(x,snr);%加入高斯白噪声
R =x*x'/N;%自相关函数
Rr=real(R);%实部
Ri=imag(R);%虚部
N=M;
J=zeros(N,N);
for k=1:N
    J(k,N+1-k)=1;
end
J=J;
T=Rr+J*Ri;  %构建新的实合成矩阵
[V,D]=eig(T);  %特征值分解
Us=V(:,M-P+1:M);
U1=Us(1:M-1,:);
U2=Us(2:M,:);% 两个方阵张成的两个子空间
[p,q]=eig(pinv(U1)*U2); % 利用最小二乘法 求得旋转不变关系矩阵，然后进行特征分解
for i=1:P;
    alpha(i)=real(asin(j*(log(q(i,i)))*lambda/(-2*pi*d))*180/pi);
end
alpha
stem(alpha,ones(1,P),'filled');grid;
axis([-90 90 0 2]);
text(alpha(3)-4,1.1,num2str(alpha(3)));text(alpha(3)-5,1.3,'信号1');
text(alpha(2)-4,1.1,num2str(alpha(2)));text(alpha(2)-5,1.3,'信号2');
text(alpha(1)-4,1.1,num2str(alpha(1)));text(alpha(1)-5,1.3,'信号3');
ylabel('DOA估计的角度值  ');
xlabel('角度');
title('ESPRIT算法DOA估计');