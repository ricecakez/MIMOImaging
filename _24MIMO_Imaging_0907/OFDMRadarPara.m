function [ f0,df,N,B,fs,tb,tc,T,K,Tp,c,PRI] = OFDMRadarPara
%PARAMETERSSYSTEM Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;              %���٣���λΪ��/��
f0 = 10e9;          %�����ź����Ƶ��
df = 1e6;
N = 64;              %number of subcarriers
B = N*df;
L = 1;                %over sampling ratio
fs = L*B;               %sampling rate
tb = 1/df;            %OFDM bit duration
alpha = 0.5;          %CP ratio
tc = alpha*tb;        %CP duration
T = tb+tc;            %transmitted OFDM symbol duration
K = 128;
Tp = K*T;             %�������ʱ�䣬��λΪ��
PRI = 1e-3;
end

