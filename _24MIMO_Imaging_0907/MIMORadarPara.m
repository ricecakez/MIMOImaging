function [Fc,Bw,fs,N,K,df,tb,tc,T,PRF,Nt,Nr,dt,dr] = MIMORadarPara
%% antenna parameters (colocated ULA)
Fc = 10e9;
c = 3e8;
lambda = c/Fc;
Nt = 4;
Nr = 16;
dt = 2;
dr = Nt*dt;

%% waveform parameters
df = 1e6;
L = 1;
N0 = 64;
N = Nt*N0;
Bw = N*df;
fs = L*Bw;
K = 128;
df = Bw/N;
tb = 1/df;
alpha = 0.5;
tc = alpha*tb;
T = tb+tc;
PRF = 1/T;

