function [ L0,Omega0,dOmega,V0,a0] = ParametersTarget

L0 = 5000;        %目标坐标系中心X坐标，默认Y坐标为0，单位为米
Omega0 = 1;       %平台转动角速度，单位为弧度
dOmega = 0;      %平台转动角速度变化率
V0 = 0;             %目标初速度
a0 = 0;              %目标加速度

end
