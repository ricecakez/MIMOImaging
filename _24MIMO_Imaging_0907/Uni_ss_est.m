function [Es] = Uni_ss_est(X,R)
    P = size(X,1);
    if mod(P,2)==0
        m = P/2;
        X1 = X(1:m,:);
        X2 = X(m+1:end,:);
        X21 = flipud(X2);
        T1 = real(X1+X21);
        T2 = -imag(X1+X21);
        T3 = imag(X1-X21);
        T4 = real(X1-X21);
        Z = [T1 T2
            T3 T4];
        [U,~,~] = svd(Z);
        U1 = U(1:m,1:R);
        U2 = U(m+1:end,1:R);
        Es = [U1+1i*U2;flipud(U1)-1i*flipud(U2)]/sqrt(2);
%         Un1 = U(1:m,(1+R):end);
%         Un2 = U(m+1:end,(1+R):end);
%         En = [Un1+1i*Un2;flipud(Un1)-1i*flipud(Un2)]/sqrt(2);
    else
        m = (P-1)/2;
        X1 = X(1:m,:);
        x = X(m+1,:); 
        X2 = X(m+2:end,:);
        X21 = flipud(X2);
        T1 = real(X1+X21);
        T2 = -imag(X1+X21);
        T3 = imag(X1-X21);
        T4 = real(X1-X21);
        Z = [T1 T2
            sqrt(2)*real(x) -sqrt(2)*imag(x) 
            T3 T4];
        [U,~,~] = svd(Z);
        U1 = U(1:m,1:R);
        u = U(m+1,1:R);
        U2 = U(m+2:end,1:R);
        Es = [U1+1i*U2;sqrt(2)*u;flipud(U1)-1i*flipud(U2)]/sqrt(2);
%         Un1 = U(1:m,(1+R):end);
%         un = U(m+1,(1+R):end);
%         Un2 = U(m+2:end,(1+R):end);
%         En = [Un1+1i*Un2;sqrt(2)*un;flipud(Un1)-1i*flipud(Un2)]/sqrt(2);
    end   
end