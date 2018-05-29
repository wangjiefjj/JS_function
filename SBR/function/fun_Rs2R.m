function [ R ] = fun_Rs2R( H, Rs )
%FUN_RS2R �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ָ���߶Ⱥ�б�����Ӧ�ĵؾ�
Re = 6373; %flat earth radius��km
t1 = Re^2+(Re+H)^2 - Rs.^2;
t2 = 2*Re*(Re+H);
R = acos(t1./t2).*Re;

end

