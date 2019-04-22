function J = permuMat(N1,N2)
    J = zeros(N1*N2);
    for n1 = 1:N1
        for n2 = 1:N2
            J((n1-1)*N2+n2,(n2-1)*N1+n1) = 1;
        end
    end
end