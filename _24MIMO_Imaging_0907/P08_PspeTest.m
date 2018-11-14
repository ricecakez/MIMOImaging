clear;clc;close all;
n = 0:99;   
s = exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);
X = corrmtx(s,12,'modified'); 
peig(X,3,'whole')
figure
pmusic(X,3,'whole')
figure
pwelch(s);
% plot(10*log10(pxx))

rng default

fs = 1000;
t = 0:1/fs:5-1/fs;

noisevar = 1/4;
x = cos(2*pi*100*t)+sqrt(noisevar)*randn(size(t));

[pxx,f] = pwelch(x,500,300,500,fs,'centered','power');

plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid

A = [1 -2.7607 3.8106 -2.6535 0.9238];
[H,F] = freqz(1,A,[],1);
plot(F,20*log10(abs(H)))

xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

rng default

x = randn(1000,1);
y = filter(1,A,x);
[Pxx,F1] = pburg(y,4,1024,1);

hold on
plot(F1,10*log10(Pxx))
legend('True Power Spectral Density','pburg PSD Estimate')