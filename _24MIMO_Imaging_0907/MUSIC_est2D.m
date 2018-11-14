function s = MUSIC_est2D(x,p,t,f,p1,p2)
%x:input signal
%p:number of input signals
%t:fast time indexes
%f:frequency indexes
%p1,p2:win length of two dimension

%% spatial smoothing
[N,K] = size(x);
nwin = p1*p2;
m = 1;
Rxx = zeros(nwin,nwin);
J = fliplr(eye(nwin));
for n = 1:N-p1+1
    for k = 1:K-p2+1
        xm = reshape(x(n:n+p1-1,k:k+p2-1),[nwin,1]);
        Rm = xm*xm';
        Rxx = Rxx + Rm;% + J*conj(Rm)*J;
        m = m+1;
    end
end
Rxx = Rxx/(m-1);
% xin = buffer(x_c,nwin,nwin-1,'nodelay');
% xin = xin'./sqrt(Lx-nwin);
[eigenvects,S] = eig(Rxx);
[eigenvals,ind]=sort(diag(S),'descend');

thresh = p(2)*eigenvals(end);
indx = find(eigenvals > thresh);
if ~isempty(indx)
    p_eff = min( p(1), length(indx) );
else
    p_eff = p(1);
end

signal_eigenvects = eigenvects(:,ind(1:p_eff));
noise_eigenvects = eigenvects(:,ind(p_eff+1:end));

s = zeros(length(t),length(f));
for i = 1:length(t)
    for j = 1:length(f)
        at = exp(-1i*2*pi*t(i)*(0:nwin-1));
        af = exp(-1i*2*pi*f(i)*(0:nwin-1));
        a = kr(at,af);
        s(i,j) = abs((a*a')/(a*noise_eigenvects*noise_eigenvects'*a'));
    end
end
end