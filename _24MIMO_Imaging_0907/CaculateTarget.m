function targetCondition = CaculateTarget( target,N_target,L0,Omega,t )
%CACULATETARGET Summary of this function goes here
%   Detailed explanation goes here
%��Ŀ����м��㣬�õ��������һ�л��Ǻ���ɢ��ϵ�����ڶ������״�����ϵ��
%���״�ľ��룬������������ת����ɵ����״����߷�����ٶ�
targetPosition = zeros(N_target,2); %����Ŀ�����ʱ��tʱ��λ��

targetPosition(:,1) = target(:,2)*cos(Omega*t)-target(:,3)*sin(Omega*t); %u=x*cos(wt)-y*sin(wt)
targetPosition(:,2) = target(:,3)*cos(Omega*t)+target(:,2)*sin(Omega*t); %v=x*sin(wt)+y*cos(wt)

targetCondition = zeros(N_target,2);

targetCondition(:,1) = target(:,1);  %Ŀ���ĺ���ɢ��ϵ��
targetCondition(:,2) = ((targetPosition(:,1)).^2+((L0+targetPosition(:,2)).^2)).^0.5; %Ŀ������״�ľ���


end

