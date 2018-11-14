function ms = sort_matrix(m,method,dim)
ms = zeros(size(m));
if dim == 1
    N = size(m,2); 
    [ms(:,1),inds] = sort(m(:,1),method);
    for n = 2:N
        ms(:,n) = m(inds,n);
    end
else
    M = size(m,1);
    [ms(1,:),inds] = sort(m(1,:),method);
    for n = 2:M
        ms(n,:) = m(n,inds);
    end
end