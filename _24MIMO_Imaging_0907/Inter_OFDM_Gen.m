function Inter_OFDM_Gen(Nt)
%% interleaved OFDM
B = 500e6;
L = 1;
fs = L*B;
N = 512;
df = B/N;
tb = 1/df;
alpha = 0.5;
tc = alpha*tb;
T = tb+tc;
K = 64;
N0 = N/Nt;   %每一个发射信号包含载频数
a = exp(1i*2*pi*ceil(4*rand(N,K))/4);
wav = zeros(N*L*(1+alpha),K,Nt);
t = 0:1/fs:(T-1/fs);
t = t(:);
eps0 = 1e-8;
for k = 1:K
    for n = 1:N0
        for nt = 1:Nt
            wav(:,k,nt) = wav(:,k,nt) + a((n-1)*Nt+nt,k)*exp(1i*2*pi*((n-1)*Nt+nt-1)*df*(t-tc));
        end
    end
end
save('InterOFDM_wav.mat');