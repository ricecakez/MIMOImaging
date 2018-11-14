function Mat = SampMat(M,N,d,method)
Mat = [];
if method == 1
    for n = 1:N
        Mat = blkdiag(Mat,[eye(M-d) zeros(M-d,d)]);
    end
else
    for n = 1:N
        Mat = blkdiag(Mat,[zeros(M-d,d) eye(M-d)]);
    end
end