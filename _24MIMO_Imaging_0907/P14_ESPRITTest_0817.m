clc;clear all,close all
N=1000;%样本数
% load un.mat;
noise=(rand(1,N)+1i*rand(1,N))/sqrt(2);
n=1:1:N;
un(n)=exp(1i*0.1*pi*n)+exp(-1i*0.3*pi*n)+noise;
%% 构造自相关矩阵
M=8;
for k=1:N-M+1
xs(:,k)=un(k+M-1:-1:k).'; %从 k+m-1 一直减到 k 并变为列向量
end
K = 2;
f = ESPRIT(un,M,K)
Rxx=xs(:,1:end-1)*xs(:,1:end-1)'/(N-M+1);
Rxy=xs(:,1:end-1)*xs(:,2:end)'/(N-M+1);
%% 利用奇异分解找到最小特征值
[V,D]=svd(Rxx);
ev=diag(D);emin=ev(end);
%% 构造矩阵对,并对其进行特征分解
z=[zeros(M-1,1),eye(M-1);0,zeros(1,M-1)];
Cxx=Rxx-emin*eye(M);
Cxy=Rxy-emin*z;
[V,D]=eig(Cxx,Rxy);
z=diag(D);
w=angle(z)/(2*pi);
err=abs(abs(z)-1);
[Enew,ad] = sort(err); %对其排序，小到大排列
w(ad([1:2]))