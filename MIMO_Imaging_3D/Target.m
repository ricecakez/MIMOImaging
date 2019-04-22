function [target,N_target] = Target
fid=fopen('Target_1.txt');
target = zeros(1,4);
i = 0;
while ~feof(fid)
    i = i+1;
    [temp_target,Count]=fscanf(fid,'%f',4);
    if Count > 0
        target(i,:) = temp_target;
    end
end
N_target = size(target,1);
fclose(fid);
end
