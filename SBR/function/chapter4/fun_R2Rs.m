function [ Rs ] = fun_R2Rs( H,R )
%FUN_RS2R �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ָ���߶Ⱥ͵ؾ����Ӧ��б��
Re = 6373e3; %flat earth radius��m
Rs = sqrt(Re.^2 + (Re+H).^2 - 2*Re.*(Re+H).*cos(R./Re));
end

