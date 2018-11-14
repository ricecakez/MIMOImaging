function [signal_return,a ]=ReturnSimulate
%RETURNSIMULATE Summary of this function goes here
%   Detailed explanation goes here
%本函数用来进行回波的模拟
%%
%获得原始数据及相关数据计算
[ L0,Omega0,dOmega,V0,a0] = ParametersTarget();
[ F0,F_sample,B,N,df,tb,tc,T,PRF,K,Tp,c] = ParametersSystem();
[target,N_target] = Target();
a = exp(1i*2*pi*ceil(4*rand(N,K))/4);
T_sample = 1/F_sample;  %系统采样时间，单位为秒
L_range = c*tc/2;
% L_range = 20;
L_min = L0-L_range/2;     %最早接受数据对应的目标距离，也即最小的接收数据的距离
L_max = L0+L_range/2;     %最晚接收数据对应的目标距离，也即最大的接收数据的距离
T_min = 2*L_min/c;      %最早的接收数据时间
T_max = 2*L_max/c;      %最晚的接收数据时间
i_symlength = round(T/T_sample); %发射OFDM符号序列长度
% i_pulse = 1:i_pulselength; %发射脉冲的脉冲序列
% t_sym = T_min + (0:i_symlength-1).'*T_sample; %发射脉冲的时间采样信号
i_rec = round((T_max-T_min)/T_sample)+i_symlength; %每次脉冲的接收窗采样序列长度
t_rec = T_min+(0:i_rec-1).'*T_sample;
%%
%初始化
signal_return = zeros(i_rec,K);

%生成回波数据
eps = 1e-8;
h1 = waitbar(0,'生成数据');
for k = 1:K
    L = L0+V0*((k-1)/PRF)+0.5*a0*((k-1)/PRF)^2;            %计算目标位置
    Omega =Omega0+dOmega*((k-1)/PRF);
    targetCondition = CaculateTarget(target,N_target,L,Omega,(k-1)/PRF);    
    for i = 1:N_target
        tau = 2*targetCondition(i,2)/c;                                      %这里是从T_min开始也就是从接收回波开始计时后的时间，为了与回波接收的数据结相吻合
%         i_return = round((t_return-T_min)/T_sample);                                              %这里计算的是在第几个采样脉冲时开始接收到回波信号
        signal_return_target = zeros(i_rec,1);
        for n = 1:N
            signal_return_target = signal_return_target + a(n,k)*exp(1i*2*pi*F0*(t_rec-tau)).*exp(1i*2*pi*(n-1)*df*(t_rec-tc-tau))...
                .*rectpuls((t_rec-tau)/T-1/2+eps);                  %这是每个目标的回波
        end
        signal_return(:,k) = signal_return(:,k)+targetCondition(i,1)*signal_return_target; %把每个目标的回波叠加进脉冲的接收回波
    end   
    waitbar(k/K);
end
delete(h1);
signal_return = awgn(signal_return,5,'measured');
% signal_return = complex(signal_return_real,signal_return_imag);
%%
%显示结果
%{
figure
plot(real(signal_return(1,:)));
title('第1次脉冲的回波信号');
figure
plot(real(signal_reference));
title('参考信号');
%}
%%
%将结果进行格式匹配
save('ReturnSignal.mat');
end

