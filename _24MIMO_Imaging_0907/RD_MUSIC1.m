function [t_e,f_e] = RD_MUSIC1(r,p,t,f)
%% inputs
%r:corrlation matrix
%p:number of input signals
%t:time indexes
%f:frequency indexes
[V,D] = eig(r);
[ds,ind] = sort(diag(D),'descend');
Vs = V(:,ind);

thresh = ds(end)*p(2);
indx = find(ds > thresh);
if ~isempty(indx)
    p_eff = min( p(1), length(indx) );
else
    p_eff = p(1);
end
En = Vs(:,p_eff+1:end);
N = length(t);
K = length(f);
L = size(En,1);
% s = zeros(1,p_eff);
pp = 0;
for k = 1:K
    af = exp(1i*2*pi*f(k)*(0:L-1));
    af = af(:);
    Q = kron(af,eye(N))'*En*En'*kron(af,eye(N));
    Q1 = inv(Q);
    if pp < p_eff 
        pp = pp + 1;
        f_e(pp) = f(k);
        s(pp) = Q1(1,1);
        at(:,pp) = Q1(:,1)/Q1(1,1);
        [s,ind] = sort(s,'descend');
        f_e = f_e(ind);
        at = at(:,ind);
    else
        if Q1(1,1) > s(end)
            s(end) = Q1(1,1);
            at(:,end) = Q1(:,1)/Q1(1,1);
            f_e(end) = f(k);
            [s,ind] = sort(s,'descend');
            f_e = f_e(ind);
            at = at(:,ind);
        end
    end
end
P = [ones(N,1) 2*pi*(0:(N-1)).'];
P1 = inv(P.'*P)*P.';
for pp = 1:p_eff
    g = -angle(at(:,pp));
    c = P1*g;
    t_e(pp) = c(2);
end
