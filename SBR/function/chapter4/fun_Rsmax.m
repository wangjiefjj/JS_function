function [ Rsmax ] = fun_Rsmax( H )
%FUN_RSMAX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%�߶Ⱥ�����б��
Re = 6373e3; %flat earth radius��m
Rsmax = sqrt(H.*(2*Re+H));%��ʽ��4.13��

end

