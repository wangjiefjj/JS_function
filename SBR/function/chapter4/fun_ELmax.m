function [ ELmax ] = fun_ELmax( H )
%FUN_ELMAX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%����H�����ĸ�����
%%
Re = 6373; %flat earth radius��km
ELmax = pi/2 - (acos(1./(1+H/Re))); %��ʽ��4.14��
ELmax = ELmax/pi*180;
end

