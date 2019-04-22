function [targetRCS,targetPos,targetVel,targetIniPos] = TargetPara
[target, N_target] = Target;
targetIniPos = [1000 1200 1500];
targetRCS = target(:,1);
targetPos = [target(:,2:4)] + ones(N_target,1)*targetIniPos;
targetVel = zeros(size(targetPos));