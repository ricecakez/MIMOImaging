clear;clc;close all;

fs = 8.0e3;
t = (0:1/fs:1).';
x1 = cos(2*pi*t*300);
x2 = cos(2*pi*t*400);
array = phased.ULA('NumElements',10,'ElementSpacing',1);
array.Element.FrequencyRange = [100e6 300e6];
fc = 150e6;
x = collectPlaneWave(array,[x1 x2],[10 20;45 60]',fc);
noise = 0.1/sqrt(2)*(randn(size(x)) + 1i*randn(size(x)));

estimator = phased.ESPRITEstimator('SensorArray',array,...
    'OperatingFrequency',fc);
doas = estimator(x + noise);
az = broadside2az(sort(doas),[20 60])