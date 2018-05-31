function [ EL, R, Rs ] = fun_ELAngle( H,R)
%FUN_ELANGLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ��������ؾ�Ĺ�ϵ����ʽ��4.7��
%%
Re = 6373; %flat earth radius��km
if nargin==1
    [R,Rs]=fun_RsR(H);
    EL=(asin(Re./Rs.*sin(R./Re)));
    EL = EL/pi*180;
else
    Rs = fun_R2Rs(H,R);
    EL=(asin(Re./Rs.*sin(R./Re)));
    EL = EL/pi*180;
end


end

