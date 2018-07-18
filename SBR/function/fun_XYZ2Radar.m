function [ Xr,Yr,Zr ] = fun_XYZ2Radar( alpha1,beta1,eta,h,Xe,Ye,Ze)
%FUN_XYZ2RADAR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%��������ϵת�״�����ϵ
% alpha1: SBR�ľ���deg,����Ϊ��������Ϊ��
% beta1�� SBR��γ��deg ����γΪ������γΪ��
% eta�� SBR�Ĺ�����deg
% h �״�߶ȣ�Re+H��
% Xe:�״��ڵع�����ϵ������
alpha1 = alpha1/180*pi;
beta1 = beta1/180*pi;
eta = eta/180*pi;
mu = asin(sin(beta1)/sin(eta));
Ar2e = [-sin(mu),         0,       -cos(alpha1)*cos(beta1);
        cos(eta)*cos(mu), sin(beta1)*sin(mu)+cos(alpha1)*cos(beta1)*sin(eta)*cos(mu),  -sin(alpha1)*cos(beta1);
        sin(eta)*cos(mu), -(sin(alpha1)*sin(mu)+cos(alpha1)*cos(eta)*cos(mu)), - sin(beta1)];
Ae2r = Ar2e.';
Z = Ae2r * [Xe;Ye;Ze-h];
Xr = Z(1);
Yr = Z(2);
Zr = Z(3);
end
