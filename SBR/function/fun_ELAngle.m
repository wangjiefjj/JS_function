function [ EL, R, Rs ] = fun_ELAngle( H )
%FUN_ELANGLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ��������ؾ�Ĺ�ϵ����ʽ��4.7��
%%
Re = 6373; %flat earth radius��km
[R,Rs]=fun_RsR(H);
EL=abs(asin(Re./Rs.*sin(R./Re)));
EL = EL/pi*180;
end

