function [ crabA ] = fun_CrabAngle( alpha1,eta, H )
%FUN_CRAB �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ƫ���Ǽ���
%% ����˵��
% alpha1: �״���Ϣ���γ�� deg
% eta���״��������н� deg
% H: �״�������µ�߶� m
%% ����delta 4.61
Re = 6373e3; %flat earth radius��m
Ve = 0.4651e3; %����ٶȣ�m/s
Vp = fun_Vp(H); %m/s
delta = Ve./Vp.*(1+H./Re);
%% ����ƫ���� 4.72
crabA = atan(delta.*sqrt(cos(alpha1/180*pi).^2-cos(eta/180*pi).^2)./(1-delta.*cos(eta/180*pi)));
end

