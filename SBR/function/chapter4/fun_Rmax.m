function [ Rmax ] = fun_Rmax( H )
%FUN_RMAX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%�߶Ⱥ����ĵؾ�
Re = 6373e3; %flat earth radius��m
Rmax = Re*acos(1./(1+H./Re));%��ʽ��4.12��
end

