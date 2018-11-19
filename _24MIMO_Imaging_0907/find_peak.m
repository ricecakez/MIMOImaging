function [inds,f] = find_peak(f,K)
for k = 1:K
    inds(k) = find(f == max(f),1);
    m = inds(k);
    n = inds(k);
    while (m>1)&&(f(m-1) <= f(m))
        m = m-1;
    end
    while (n<length(f))&&(f(n+1) <= f(n))  
        n = n+1;
    end
    f(m:n) = 0;
end