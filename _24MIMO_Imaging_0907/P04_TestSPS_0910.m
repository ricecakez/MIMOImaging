clear;clc;close all;

N = 10;
d = 0.5;
elementPos = (0:N-1)*d;
angles = [0 -25];
ac = [1 1/5];
scov = ac'*ac;
R = sensorcov(elementPos,angles,db2pow(-5),scov);

Nsig = 2;
doa =  rootmusicdoa(R,Nsig)

Nsig = 2;
L = 2;
RSM = spsmooth(R,L);
doasm = rootmusicdoa(RSM,Nsig)