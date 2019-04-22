function [f1,f2,f3] = MDF_3D(X,R)
    I = size(X);
    for n = 1:3
        if mod(I(n),2) == 0
            I1(n) = I(n)/2;
            I2(n) = I(n)/2+1;
        else
            I1(n) = (I(n)+1)/2;
            I2(n) = (I(n)+1)/2;
        end
    end
    for i11 = 1:I1(1)
        for i12 = 1:I2(1)
            for i21 = 1:I1(2)
                for i22 = 1:I2(2)
                    for i31 = 1:I1(3)
                        for i32 = 1:I2(3)
                            u = ((i11-1)*I1(2)+i21-1)*I1(3)+i31;
                            v = ((i12-1)*I2(2)+i22-1)*I2(3)+i32;
                            X1(u,v) = X(i11+i12-1,i21+i22-1,i31+i32-1);
                            Y1(u,v) = conj(X(I(1)-(i11+i12-1)+1,I(2)-(i21+i22-1)+1,I(3)-(i31+i32-1)+1));
                        end
                    end
                end
            end
        end
    end
    X_ = [X1;Y1];
    [U,~,~] = svds(X_,R);
    M = size(U,1);
    tmp = U(1:M/2,:)\U((M/2+1):end,:);
    [T_1,~]= eig(tmp);
    G = U(1:M/2,:)*T_1;
    J1 = kron([eye(I1(1)-1) zeros(I1(1)-1,1)],eye(I1(2)*I1(3)));
    J2 = kron([zeros(I1(1)-1,1) eye(I1(1)-1)],eye(I1(2)*I1(3)));       
    J3 = kron(eye(I1(1)),kron([eye(I1(2)-1) zeros(I1(2)-1,1)],eye(I1(3))));
    J4 = kron(eye(I1(1)),kron([zeros(I1(2)-1,1) eye(I1(2)-1)],eye(I1(3))));
    J5 = kron(eye(I1(1)*I1(2)),[eye(I1(3)-1) zeros(I1(3)-1,1)]);
    J6 = kron(eye(I1(1)*I1(2)),[zeros(I1(3)-1,1) eye(I1(3)-1)]);      
    f1 = angle(diag((J1*G)\(J2*G)))/2/pi;
    f2 = angle(diag((J3*G)\(J4*G)))/2/pi;
    f3 = angle(diag((J5*G)\(J6*G)))/2/pi;
    