function [ F0,F_sample,B,N,df,tb,tc,T,PRF,K,Tp,c] = ParametersSystem
%PARAMETERSSYSTEM Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;              %光速，单位为米/秒
F0 = 10e9;          %发射信号最低频率
% F_sample = 5e9;     %系统采样频率（这里F_sample需<Fc，否则在Keystone函数的调用中会出错）
% B = 100e6;            %信号带宽
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
Tp = K*T;             %脉冲持续时间，单位为秒


end

