function [targetRCS,targetPos,targetVel,targetIniPos] = TargetPara
[target, N_target] = Target;
targetIniPos = [3000 4000 1000];
targetRCS = target(:,1);
targetPos = [target(:,2:3) zeros(N_target,1)] + ones(N_target,1)*targetIniPos;
targetVel = zeros(size(targetPos));