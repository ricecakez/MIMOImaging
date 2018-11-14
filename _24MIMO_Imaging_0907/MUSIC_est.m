function s = MUSIC_est(x,p,f,nwin)
%x:input signal
%p:number of input signals
%t:threshhold of the eigenvalue
%f:frequency indexes
%nwin:win length

%% spatial smoothing
Lx = length(x);
xin = buffer(x,nwin,nwin-1,'nodelay');
xin = xin'./sqrt(Lx-nwin);
[~,S,eigenvects] = svd(xin,0);
eigenvals = diag(S).^2; 

thresh = p(2)*eigenvals(end);
indx = find(eigenvals > thresh);
if ~isempty(indx)
    p_eff = min( p(1), length(indx) );
else
    p_eff = p(1);
end

signal_eigenvects = eigenvects(:,1:p_eff);
noise_eigenvects = eigenvects(:,p_eff+1:end);

s = zeros(size(f));
for i = 1:length(f)
    a = exp(-1i*2*pi*f(i)*(0:nwin-1));
    s(i) = abs((a*a')/(a*noise_eigenvects*noise_eigenvects'*a'));
end
end