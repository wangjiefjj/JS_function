function [ GA,R,Rs ] = fun_GrazeAngle( H, R, Rs)
%FUN_ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �ӵؽ���ؾ�Ĺ�ϵ����ʽ��4.6��
% H:����״�����µ�ľ���,m
%%
Re = 6373e3; %flat earth radius��m
if nargin==1
    [R,Rs] = fun_RsR(H);
end
GA = (acos((Re+H)./Rs .* sin(R/Re)));
GA = rad2deg(GA);
end

