function [ Vp ] = fun_Vp( H )
%FUN_VP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%H���߶ȵ�λΪm
%% �������Ǹ߶ȼ��������ٶ�
G = 6.673e-11; %%���泣��
Me = 5.965e24; %%��������kg
Re = 6373*10^3; %flat earth radius��m
Vp = sqrt(G*Me./(Re+H));%m/s
end

