function [ FootLength ] = fun_FootLength( H,EL_beam )
%FUN_FOOTLENGTH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% EL_beam: ����������ȡ�rad
% H:����״�����µ�ľ���,km
%������ӡ����
Re = 6373e3; %flat earth radius��m
GrazeT = fun_GrazeTH(H,EL_beam,1);
GrazeH = fun_GrazeTH(H,EL_beam,2);
FootLength = Re.*(GrazeH-GrazeT-EL_beam);

end

