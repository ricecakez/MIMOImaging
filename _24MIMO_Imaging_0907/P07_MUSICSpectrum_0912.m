clear;clc;close all;

n = 0:63;
% L = length(n);
d = 0.1e-2;
L = round(1/d);
s = exp(1i*2*pi*d*n) + exp(1i*2*pi*5*d*n) + 0.1*randn(size(n));%+ exp(1i*pi/3*n) + randn(1,100);
figure
plot(linspace(0,1,L),abs(fft(s,L)))
xlim([0,1e-1])
% Hs = spectrum.music(2,40);
% figure
% pseudospectrum(Hs,s)
% 
f = linspace(0,1e-2,2048);
[S,w] = pmusic(s,[2,40],f,[],40);
figure
plot(w,mag2db(S))
grid on

S = MUSIC_est(s,[2,5],f,40);
figure
plot(f,mag2db(S))
