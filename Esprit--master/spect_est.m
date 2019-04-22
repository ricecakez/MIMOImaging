function [f] = spect_est(X,R,iterNum,thresh)
U = CP_ALS(X,R,iterNum,thresh);
N = length(U);
for n = 1:N
    Mn = size(U{n},1);
    Pn = [ones(Mn,1) (0:Mn-1).'];
    Pn1 = pinv(Pn);
    for r = 1:R
        hn = phase(U{n}(:,r))/2/pi;
        cn = Pn1*hn;
        f(n,r) = cn(2);
    end
end