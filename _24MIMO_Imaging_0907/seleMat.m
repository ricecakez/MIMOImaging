function [K1,K2] = seleMat(M)
    m = M-1;
%     J1 = [eye(m),zeros(m,1)];
    J2 = [zeros(m,1),eye(m)];
    Qm = UniMat(m);
    QM = UniMat(M);
    tmp = Qm'*J2*QM;
    K1 = real(tmp);
    K2 = imag(tmp);
end