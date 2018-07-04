function [ R,Rs ] = fun_RsR( H )
%FUN_RSR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵����
%% Э��Rs�͵ؾ�R�Ĺ�ϵ����ʽ��4.3��
% H:����״�����µ�ľ���,km
%%
Re = 6373e3; %flat earth radius��m
R = 0:1e3:fun_Rmax(H);
Rs = sqrt(Re^2 + (Re+H)^2 - 2*Re*(Re+H)*cos(R/Re));


end

