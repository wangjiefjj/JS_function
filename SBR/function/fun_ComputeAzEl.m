function [All] = fun_ComputeAzEl( alpha1,beta1,h1,eta,alpha2,beta2,h2 )
%FUN_COMPUTEAZEL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%  ����Ŀ�꾭γ�ߺ��״ﾭγ�߹����������㲨����ָ��
%alpha1: �״ﾭ��
%beta1: �״�γ��
%h1: �״�ߣ�H��
%eta: �״������
Re = 6373e3;
[X1,Y1,Z1] = fun_JWH2XYZ(alpha1,beta1, h1+Re); %���Ǿ�γ��תXYZ
[X2,Y2,Z2] = fun_JWH2XYZ(alpha2,beta2, h2+Re);   %Ŀ��㾭γ��תXYZ
All.Rs = sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);
All.R = fun_Rs2R(h1,All.Rs);
%���㲨��ָ��ĸ�����
All.graze = fun_GrazeAngle(h1,All.R,All.Rs);
All.El = fun_ELAngle(h1,All.R);
% �����״�����ϵ����
All.mu = asin(sin(beta1/180*pi)/sin(eta/180*pi));
[All.Xr1,All.Yr1,All.Zr1] = fun_XYZ2Radar(alpha1,beta1,eta,Re+h1, X1,Y1,Z1);
% Xr1 = -sin(mu)*X1+cos(eta/180*pi)*cos(mu)*Y1+sin(eta/180*pi)*cos(mu)*(Z1-Re-H);
All.Az = abs(acos(All.Xr1/All.Rs/sin(All.El/180*pi)))/pi*180;
end

