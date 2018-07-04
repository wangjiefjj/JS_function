function [ GrazeTH] = fun_GrazeTH(H, EL_beam,opt )
%FUN_ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% EL_beam: ����������ȡ�rad
%% ����Զ���ص������
Re = 6373e3; %flat earth radius��m
[EL,R,Rs] = fun_ELAngle(H); 
EL = EL/180*pi;
if opt==1 %%Զ�ص�
    temp = (1+H/Re)*sin(EL+EL_beam/2);
    GrazeTH = abs(acos(temp));
else
	temp = (1+H/Re)*sin(EL-EL_beam/2);
	GrazeTH = abs(acos(temp));  
end


end

