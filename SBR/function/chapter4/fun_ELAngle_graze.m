function [EL] = fun_ELAngle_graze(H,graze)
%FUN_ELANGLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ��������ؾ�Ĺ�ϵ����ʽ��4.8��
%% H:�߶�m��graze���ӵؽ�deg,EL:������deg
Re = 6373e3; %flat earth radius��m
EL = asin(1/(1+H/Re)*cos(graze/180*pi))/pi*180;

end

