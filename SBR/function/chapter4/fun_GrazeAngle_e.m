function [ graze ] = fun_GrazeAngle_e( H, R, Rs,alpha2)
%FUN_ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �ӵؽ���ؾ�Ĺ�ϵ��Բ���򣬹�ʽ��4B.12��(4B.34)
% H:����״�����µ�ľ���,m
% R:�ؾ�m
% Rs:б��m
% alpha2��Ŀ��γ��deg
%%
e = 0.08199;%���������� 
graze = fun_GrazeAngle(H,R,Rs); %Բ���������deg
v = atan(e^2*sin(alpha2)/(2*(1-e^2*cos(alpha2)^2)))/pi*180; %������deg
graze = graze + v;
end

