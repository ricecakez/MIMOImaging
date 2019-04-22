function [mu,nu] = StepUESPRIT(Y,M,N,Nt,I)
[P,L] = size(Y);
if mod(P,2) == 0
    m = P/2;
else
    m = (P-1)/2;
    y = Y(m+1,:);
end
Y1 = Y(1:m,:);
Y2 = Y((end-m+1):end,:);
Y21 = IEM(m) * Y2;
if mod(P,2) == 0
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        imag(Y1-Y21) real(Y1-Y21)];
else
    Z = [real(Y1+Y21) -imag(Y1+Y21);
        sqrt(2)*real(y) -sqrt(2)*image(y);
        imag(Y1-Y21) real(Y1-Y21)];
end
[U,D,~] = svd(Z);
Es = UniMat(P)*U(:,1:I);
J1 = kron(eye(N),[eye(M-1),zeros(M-1,1)]);
J2 = kron(eye(N),[zeros(M-1,1), eye(M-1)]);
Lambda_mu = pinv(J1*Es)*(J2*Es);
[V,D] = eig(Lambda_mu);
mu = angle(diag(D))/2/pi;
% Te = rowPermutateMat(M,N);
% Es2 = Te*Es;
A = Es*V;
Nr = M/Nt;
N_c = kron(kron(ones(Nr,1),exp(-1i*2*pi*(0:Nt-1).'*mu.'/Nt)),ones(M,1));
% A = N_c.*A;
Phi = (A(1:end-N,:))\(A(N+1:end,:));
nu = angle(diag(Phi))/2/pi;

