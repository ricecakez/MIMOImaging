clear;clc;close all

%% radar parameters
f0 = 10e9;
N = 128;
B = 128e6;
df = B/N;
tb = 1/df;
K = 1000;
c = 3e8;
L = 1;
fs = B;
dR = c/2/fs;
alpha = 0.5;
tc = alpha*tb;
T = tb+tc;
Tp = T*K;
lambda = c/f0;
dv = lambda/2/Tp;
Ruam = T*c/2;
Ns = N*L*(1+alpha);
Nc = N*L*alpha;
wav = zeros(Ns,K);
a = exp(1i*2*pi*ceil(rand(N,K)*4)/4);   %QPSK
for k = 1:K
    wav(Nc+1:end,k) = L*N*ifft(a(:,k)/sqrt(N),L*N);
    wav(1:Nc,k) = wav(N*L+1:end,k);
end

PRF = 1/T;

antenna = phased.IsotropicAntennaElement('FrequencyRange',[2e9 15e9]);      %antenna setting
radar_motion = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
R0 = 3e3;
tgtvel_xy = 500;
omega = tgtvel_xy/R0;
% R0 = 20;
theta0 = 0.78;
tgtpos_xy = [-10 -6 -6 4 4; 0 1 -1 -1 1];


D = size(tgtpos_xy,2);
tgtpos_uv = [tgtpos_xy(1,:)*cos(theta0)-tgtpos_xy(2,:)*sin(theta0);
    tgtpos_xy(1,:)*sin(theta0)+tgtpos_xy(2,:)*cos(theta0)];
tgtpos_rad = [tgtpos_uv(1,:);tgtpos_uv(2,:)+R0;zeros(1,D)];
tgtvel_rad = [tgtvel_xy*cos(theta0);tgtvel_xy*sin(theta0);0]*ones(1,D);
tgtrcs = ones(1,D);
targets = phased.RadarTarget('MeanRCS',tgtrcs,'PropagationSpeed',c,'OperatingFrequency',f0);
targets_motion = phased.Platform('InitialPosition',tgtpos_rad,'Velocity',tgtvel_rad);

transmitter = phased.Transmitter('PeakPower',50e3,'Gain',36);
receiver = phased.ReceiverPreamp('Gain',20,'NoiseFigure',5,...
    'ReferenceTemperature',290,'SampleRate',fs,...
    'EnableInputPort',true,'SeedSource','Property','Seed',1e3);
txradiator = phased.Radiator('Sensor',antenna,'OperatingFrequency',f0,...
    'PropagationSpeed',c,'WeightsInputPort',true);
rxcollector = phased.Collector('Sensor',antenna,'OperatingFrequency',f0,...
    'PropagationSpeed',c);
channel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',f0,'SampleRate',fs,'TwoWayPropagation',true);

xr = complex(zeros(Ns,K));

for k = 1:K
    % Update radar and target positions
    [radar_pos,radar_vel] = radar_motion(T);
    radar_pos = [0 2950 0].';
    [tgt_pos,tgt_vel] = targets_motion(T);
    [tgt_rang,tgt_ang] = rangeangle(tgt_pos,radar_pos);
%     tgt_rang = tgt_rang - 36e3;

    % Transmit FMCW waveform
    sig = wav(:,k);
    txsig = transmitter(sig);
    
    txsig = txradiator(txsig,tgt_ang,1);

    % Propagate the signal and reflect off the target
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    txsig = targets(txsig);

    % the received radar return
    rxsig = rxcollector(txsig,tgt_ang);
    xr(:,k) = receiver(rxsig,ones(Ns,1));   
    y(:,k) = xr(Nc+1:end,k); %ȥCP
    tmp = fft(y(:,k));  %decoding
    I(:,k) = tmp(1:N);
    I1(:,k) = I(:,k)./a(:,k);
    Y1(:,k) = ifft(I1(:,k));
    plot((0:N-1)*dR,abs(Y1(:,1)))
    hold on
end

for n = 1:N
    Y2(n,:) = ifft(Y1(n,:),1024);
end

figure
imagesc((0:1024-1)*dv*K/1024/omega,2950+(0:N-1)*dR,abs(Y2));



