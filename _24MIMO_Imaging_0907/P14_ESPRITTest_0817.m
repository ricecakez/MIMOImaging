clc;clear all,close all
N=1000;%������
% load un.mat;
noise=(rand(1,N)+1i*rand(1,N))/sqrt(2);
n=1:1:N;
un(n)=exp(1i*0.1*pi*n)+exp(-1i*0.3*pi*n)+noise;
%% ��������ؾ���
M=8;
for k=1:N-M+1
xs(:,k)=un(k+M-1:-1:k).'; %�� k+m-1 һֱ���� k ����Ϊ������
end
K = 2;
f = ESPRIT(un,M,K)
Rxx=xs(:,1:end-1)*xs(:,1:end-1)'/(N-M+1);
Rxy=xs(:,1:end-1)*xs(:,2:end)'/(N-M+1);
%% ��������ֽ��ҵ���С����ֵ
[V,D]=svd(Rxx);
ev=diag(D);emin=ev(end);
%% ��������,��������������ֽ�
z=[zeros(M-1,1),eye(M-1);0,zeros(1,M-1)];
Cxx=Rxx-emin*eye(M);
Cxy=Rxy-emin*z;
[V,D]=eig(Cxx,Rxy);
z=diag(D);
w=angle(z)/(2*pi);
err=abs(abs(z)-1);
[Enew,ad] = sort(err); %��������С��������
w(ad([1:2]))