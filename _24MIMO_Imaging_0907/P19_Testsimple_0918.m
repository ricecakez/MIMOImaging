clear;clc;close all
fs = 8e3;
fop = 3e8;
henv = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fop);
pos1 = [1000;0;0];
pos2 = [300;200;50];
vel1 = [0;0;0];
vel2 = [0;0;0];
x = 1:5;
x = x(:);
y = step(henv,x,...
    pos1,...
    pos2,...
    vel1,...
    vel2);
disp(y)

R = sqrt( (pos1-pos2)'*(pos1-pos2));
lambda = physconst('Lightspeed')/fop;
L = (4*pi*R/lambda)^2
