function [targetRCS,targetPos,targetVel,targetIniPos] = TargetParaISAR
[target, N_target] = Target;
targetIniPos = [1200;1600;0];
targetRCS = target(:,1);
targetPos = [target(:,2:3) zeros(N_target,1)].';
targetVel = [1000;0;0];