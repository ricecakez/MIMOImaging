function X1 = NormCol(X)
[Ro,Co] = size(X);
X1 = zeros(size(X));
% G = [ones(Ro,1) (0:Ro-1).'];
for cc = 1:Co
    X1(:,cc) = X(:,cc)/X(1,cc);
%     a = phase(X(:,cc))/2/pi;
%     tmp = G\a;
%     phi(cc) = tmp(2);
end
% X1 = exp(1i*2*pi*(0:Ro-1).'*phi);