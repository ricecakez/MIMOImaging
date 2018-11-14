clear;clc;close all;

Q = [0,1];
xtm = linspace(-1,1,20);
xrn = linspace(-1,1,20);
c = 3e8;
for m = 1:length(xtm)
    for n = 1:length(xrn)
        tau(m,n) = (sqrt(Q(2)^2+(xtm(m)-Q(1))^2)+sqrt((Q(1)-xrn(n))^2+Q(2)^2))/c;
    end
%     figure(1)
%     plot(xrn,tau(m,:));
%     hold on
%     axis xy
end
figure
surf(xrn,xtm,tau,'linestyle','none');

