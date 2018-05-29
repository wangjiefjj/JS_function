function [ crabA ] = fun_CrabAngle( alpha1,eta, H )
%FUN_CRAB �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ƫ���Ǽ���
%% ����˵��
% alpha1: �״���Ϣ���γ��
% eta���״��������н�
% H: �״�������µ�߶�
%% ����delta 4.61
Re = 6373; %flat earth radius��km
Ve = 0.4651; %����ٶȣ�km/s
Vp = fun_Vp(H)./1000; %km/s
delta = Ve./Vp.*(1+H./Re);
%% ����ƫ���� 4.72
crabA = atan(delta.*sqrt(cos(alpha1).^2-cos(eta).^2)./(1-delta.*cos(eta)));
end

