function [ F0,F_sample,B,N,df,tb,tc,T,PRF,K,Tp,c] = ParametersSystem
%PARAMETERSSYSTEM Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;              %���٣���λΪ��/��
F0 = 10e9;          %�����ź����Ƶ��
% F_sample = 5e9;     %ϵͳ����Ƶ�ʣ�����F_sample��<Fc��������Keystone�����ĵ����л����
% B = 100e6;            %�źŴ���
df = 1e6;

N = 64;              %number of subcarriers
B = N*df;
L = 1;                %over sampling ratio
F_sample = L*B;       %sampling rate
% df = B/N;             %subcarriers frequency interval
tb = 1/df;            %OFDM bit duration
alpha = 0.5;          %CP ratio
tc = alpha*tb;        %CP duration
T = tb+tc;            %transmitted OFDM symbol duration
K = 64;
PRF = 1/T;            %PRF
Tp = K*T;             %�������ʱ�䣬��λΪ��


end

