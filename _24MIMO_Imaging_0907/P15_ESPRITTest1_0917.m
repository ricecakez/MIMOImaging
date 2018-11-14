clear;clc;close all;

% N = 100;
% % load un.mat;
M = 24;
K = 2;
N=1000;%Ñù±¾Êý
% load un.mat;
noise=(rand(1,N)+1i*rand(1,N))/2;
n=1:1:M;
n=n(:);
un(n,:)=exp(1i*0.4*pi*n).*exp(1i*randn(1,N))+exp(1i*0.3*pi*n).*exp(1i*randn(1,N))+noise;
f = Unitary_ESPRIT_1009(un,K);