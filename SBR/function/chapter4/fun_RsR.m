function [ R,Rs ] = fun_RsR( H )
%FUN_RSR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵����
%% Э��Rs�͵ؾ�R�Ĺ�ϵ����ʽ��4.3��
% H:����״�����µ�ľ���,km
%%
Re = 6373; %flat earth radius��km
R = 0:fun_Rmax(H);
Rs = sqrt(Re^2 + (Re+H)^2 - 2*Re*(Re+H)*cos(R/Re));


end

