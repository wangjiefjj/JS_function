function [ CrabM ] = fun_CrabMagnitude( alpha1, eta, H)
%FUN_CRABMAGNITUDE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ƫ���Ǽ���
%% ����˵��
% alpha1: �״���Ϣ���γ��
% eta���״��������н�
% H: �״�������µ�߶�
%% ����delta 4.61
Re = 6378e3; %flat earth radius��m
Ve = 0.4651e3; %����ٶȣ�m/s
Vp = fun_Vp(H); %m/s
delta = Ve./Vp.*(1+H./Re);
%% ����ƫ������ 4.73
CrabM = sqrt(1+delta.^2.*cos(alpha1/180*pi).^2-2.*delta.*cos(eta/180*pi));
end

