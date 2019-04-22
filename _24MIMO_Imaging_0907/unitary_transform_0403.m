function Z = unitary_transform_0403(X)
    [M,K] = size(X);
    Z = zeros(size(X));
    if mod(M,2) == 0
        P = M/2;
        X1 = X(1:P,:);
        X2 = X(P+1:end,:);
        X21 = IEM(P)*X2;
        Z = [real(X1) + real(X21) -imag(X1)-imag(X21);
            imag(X1) - imag(X21) real(X1)-real(X21)];
    else
        P = (M+1)/2;
        X1 = X(1:P-1,:);
        x = X(P,:);
        X2 = X(P+1:end,:);
        X21 = IEM(P)*X2;
        Z = [real(X1) + real(X21) -imag(X1)-imag(X21);
            sqrt(2)*real(x) -sqrt(2)*imag(x);
            image(X1) - imag(X21) real(X1)-real(X21)];        
end
    