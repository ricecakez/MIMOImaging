function Z = UniTrans_Mat(X)
    P = size(X,1);
    if mod(P,2)==0
        m = P/2;
        X1 = X(1:m,:);
        X2 = X(m+1:end,:);
        X21 = IEM(m)*X2;
        T1 = real(X1+X21);
        T2 = -imag(X1+X21);
        T3 = imag(X1-X21);
        T4 = real(X1-X21);
        Z = [T1 T2
            T3 T4];
    else
        m = (P-1)/2;
        X1 = X(1:m,:);
        x = X(m+1,:);
        X2 = X(m+2:end,:);
        X21 = IEM(m)*X2;
        T1 = real(X1+X21);
        T2 = -imag(X1+X21);
        T3 = imag(X1-X21);
        T4 = real(X1-X21);
        Z = [T1 T2
            sqrt(2)*real(x) -sqrt(2)*imag(x) 
            T3 T4];
    end
end