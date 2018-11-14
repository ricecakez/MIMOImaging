clear;clc;close all;

M = 8;
N = 6;
K = 3;
theta = [10 20 30]*pi/180;
ar = exp(-1i*pi*(0:N-1)'*sin(theta));
Js1 = [eye(N-1) zeros(N-1,1)];
Js2 = [zeros(N-1,1) eye(N-1)];
K = N-1;
% ar1 = [Js1*fliplr(eye(N))*conj(ar(:,1));ar(:,1)];
% Js1*ar(:,1)*exp(-1i*pi*sin(theta(1))) - Js2*ar(:,1)
% Q = [eye(K) zeros(K,1) 1i*eye(K);zeros(1,K),sqrt(2),zeros(1,K);fliplr(eye(K)) zeros(K,1) -1i*fliplr(eye(K))]/sqrt(2);
% Q'*ar1
