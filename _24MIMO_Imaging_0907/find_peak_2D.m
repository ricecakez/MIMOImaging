function [m,n,Smax] = find_peak_2D(S,I)
[M,N] = size(S);
i = 1;
while i <= I
    Smax(i) = max(max(S));
    [m(i),n(i)] = find(S==Smax(i),1);
    m1 = m(i); m2 = m(i); n1 = n(i); n2 = n(i);
    while (m1>1) && (m2<M) && (n1>1) && (n2<N) 
        while (S(m1-1,n1)<=S(m1,n1))&&(S(m1-1,n2)<=S(m1,n2))&&(S(m2+1,n1)<= S(m2,n1))&&(S(m2+1,n2)<=S(m2,n2))...
                && (S(m1,n1-1)<=S(m1,n1))&&(S(m1,n2+1)<=S(m1,n2))&&(S(m2,n1-1)<=S(m2,n1))&&(S(m2+1,n2)<=S(m2,n2))...
            m1 = m1-1;
            m2 = m2+1;
            n1 = n1-1;
            n2 = n2+1;       
        end
        break;
    end
    S(m1:m2,n1:n2) = 0;
    i = i + 1;
end