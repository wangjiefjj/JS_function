function [ alpha2, beta2, h ] = fun_Radar2JWH2( Az,El,Rs, H, alpha1,beta1,phi )
%FUN_RADAR2JWH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% ֪������������Ӳ�����SBR����ϵ�е����꣬���Ӳ����ڵ����ϵľ�γ����Ϣ
% ���ݡ�����״��ʱ����Ӧ�Ӳ����Ƽ�����P23��P24����
%�͡�����״ﻷ�������붯Ŀ���⼼���о�����3.3���õ� 
%% ���� 
% Az��������λ��deg
% El������������deg
% Rs: ����Ȥ����SBR��б��m
% alpha1: SBR�ľ���deg,����Ϊ��������Ϊ��
% beta1�� SBR��γ��deg ����γΪ������γΪ��
% phi�� SBR�Ĺ�����deg
Re = 6373e3;
%% ��ʼ
Az = Az/180*pi;
El = El/180*pi;
alpha1 = alpha1/180*pi;
beta1 = beta1/180*pi;
phi = phi/180*pi;
mu = asin(sin(beta1)/sin(phi));
%%�����״�������Ϊ
xr = [1,0,0];
yr = [0,1,0];
zr = [0,0,1];
%%���Ӳ����״������е�����Ϊ
coor_c_r = Rs*(sin(El)*cos(Az)*xr + sin(El)*sin(Az)*yr + cos(El)*zr);
%%�״�ת����������Ϊ
coor_c_r = (coor_c_r+[0,0,-(Re+H)]).'
Ar2e = [-sin(mu),         0,       -cos(alpha1)*cos(beta1);
        cos(phi)*cos(mu), sin(beta1)*sin(mu)+cos(alpha1)*cos(beta1)*sin(phi)*cos(mu),  -sin(alpha1)*cos(beta1);
        sin(phi)*cos(mu), -(sin(alpha1)*sin(mu)+cos(alpha1)*cos(phi)*cos(mu)), - sin(beta1)];
%%�ڵع�����ϵ�µ�����
coor_c_e = Ar2e*coor_c_r
%%������״ﻷ�������붯Ŀ���⼼���о�����3.3���������
%%�ȼ���Ϊ��Բ����ƫ����e=0������ƫ��f=0��
e = 0;%sqrt(0.0068035111);
f = 0;%4-sqrt(4*e^2)/2;
a = 6378249.145; % ���򳤰��ᣬ�����õ�Բ����뾶m����Ϊ��ʱ��ΪΪԲ����
for i = 1:2
    p = sqrt(coor_c_e(1)^2+coor_c_e(2)^2);
    r = sqrt(p^2+coor_c_e(3)^2);
    u = atan(coor_c_e(3)/p*((1-f)+e^2*a/r));
    alpha2 = atan(coor_c_e(2)/coor_c_e(1)); %�Ӳ��ľ���
    t1 = coor_c_e(3)*(1-f)+e^2*a*sin(u)^3;
    t2 = (1-f)*(p-e^2*a*cos(u)^3);
    beta2 = atan(t1/t2); %�Ӳ���γ��
    h = p*cos(alpha2);
    v = a/sqrt(1-e^2*sin(beta2));
    coor_c_e(1) = (v+h)*cos(beta2)*cos(alpha1);
    coor_c_e(2) = (v+h)*cos(beta2)*sin(alpha1);
    coor_c_e(3) = ((1-e^2)*v+h)*sin(beta2);
end
end

