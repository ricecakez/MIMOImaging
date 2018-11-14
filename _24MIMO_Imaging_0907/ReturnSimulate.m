function [signal_return,a ]=ReturnSimulate
%RETURNSIMULATE Summary of this function goes here
%   Detailed explanation goes here
%�������������лز���ģ��
%%
%���ԭʼ���ݼ�������ݼ���
[ L0,Omega0,dOmega,V0,a0] = ParametersTarget();
[ F0,F_sample,B,N,df,tb,tc,T,PRF,K,Tp,c] = ParametersSystem();
[target,N_target] = Target();
a = exp(1i*2*pi*ceil(4*rand(N,K))/4);
T_sample = 1/F_sample;  %ϵͳ����ʱ�䣬��λΪ��
L_range = c*tc/2;
% L_range = 20;
L_min = L0-L_range/2;     %����������ݶ�Ӧ��Ŀ����룬Ҳ����С�Ľ������ݵľ���
L_max = L0+L_range/2;     %����������ݶ�Ӧ��Ŀ����룬Ҳ�����Ľ������ݵľ���
T_min = 2*L_min/c;      %����Ľ�������ʱ��
T_max = 2*L_max/c;      %����Ľ�������ʱ��
i_symlength = round(T/T_sample); %����OFDM�������г���
% i_pulse = 1:i_pulselength; %�����������������
% t_sym = T_min + (0:i_symlength-1).'*T_sample; %���������ʱ������ź�
i_rec = round((T_max-T_min)/T_sample)+i_symlength; %ÿ������Ľ��մ��������г���
t_rec = T_min+(0:i_rec-1).'*T_sample;
%%
%��ʼ��
signal_return = zeros(i_rec,K);

%���ɻز�����
eps = 1e-8;
h1 = waitbar(0,'��������');
for k = 1:K
    L = L0+V0*((k-1)/PRF)+0.5*a0*((k-1)/PRF)^2;            %����Ŀ��λ��
    Omega =Omega0+dOmega*((k-1)/PRF);
    targetCondition = CaculateTarget(target,N_target,L,Omega,(k-1)/PRF);    
    for i = 1:N_target
        tau = 2*targetCondition(i,2)/c;                                      %�����Ǵ�T_min��ʼҲ���Ǵӽ��ջز���ʼ��ʱ���ʱ�䣬Ϊ����ز����յ����ݽ����Ǻ�
%         i_return = round((t_return-T_min)/T_sample);                                              %�����������ڵڼ�����������ʱ��ʼ���յ��ز��ź�
        signal_return_target = zeros(i_rec,1);
        for n = 1:N
            signal_return_target = signal_return_target + a(n,k)*exp(1i*2*pi*F0*(t_rec-tau)).*exp(1i*2*pi*(n-1)*df*(t_rec-tc-tau))...
                .*rectpuls((t_rec-tau)/T-1/2+eps);                  %����ÿ��Ŀ��Ļز�
        end
        signal_return(:,k) = signal_return(:,k)+targetCondition(i,1)*signal_return_target; %��ÿ��Ŀ��Ļز����ӽ�����Ľ��ջز�
    end   
    waitbar(k/K);
end
delete(h1);
signal_return = awgn(signal_return,5,'measured');
% signal_return = complex(signal_return_real,signal_return_imag);
%%
%��ʾ���
%{
figure
plot(real(signal_return(1,:)));
title('��1������Ļز��ź�');
figure
plot(real(signal_reference));
title('�ο��ź�');
%}
%%
%��������и�ʽƥ��
save('ReturnSignal.mat');
end

