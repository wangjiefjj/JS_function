function [ GA,R,Rs ] = fun_GrazeAngle( H, R, Rs)
%FUN_ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �ӵؽ���ؾ�Ĺ�ϵ����ʽ��4.6��
% H:����״�����µ�ľ���,km
%%
Re = 6373; %flat earth radius��km
if nargin==1
    [R,Rs] = fun_RsR(H);
end
GA = abs(acos((Re+H)./Rs .* sin(R/Re)));
GA = GA/pi*180;
end

