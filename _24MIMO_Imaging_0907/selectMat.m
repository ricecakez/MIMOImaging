function J = selectMat(M,type)
if type == 1
    J = [eye(M-1) zeros(M-1,1)];
else
    J = [zeros(M-1,1) eye(M-1)];
end
