function [t1,r1] = sort2D(t,r,method)
if method == 1
    [t1,inds] = sort(t);
else
    [t1,inds] = sort(t,'descend');
end
r1 = r(inds);