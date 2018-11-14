function Q = UniMat(N)
if mod(N,2) == 0
    P = N/2;
    Q = [eye(P) 1i*eye(P)
        IEM(P) -1i*IEM(P)]/sqrt(2);
else
    P = (N-1)/2;
    Q = [eye(P) zeros(P,1) 1i*eye(P)
        zeros(1,P) sqrt(2) zeros(1,P)
        IEM(P) zeros(P,1) -1i*IEM(P)]/sqrt(2);
end