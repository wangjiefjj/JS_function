function [ alpha2, beta2, h] = fun_Radar2JWH( Az,El,Rs, alpha1,beta1,phi )
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
%% ��ʼ
Az = Az/180*pi;
El = El/180*pi;
alpha1 = alpha1/180*pi;
beta1 = beta1/180*pi;
phi = phi/180*pi;
mu = asin(sin(beta1)/sin(phi));
%%����ع�����ϵ��������ֱ�Ϊ
xe = [1,0,0];
ye = [0,1,0];
ze = [0,0,1];
%%���״������ϵΪ
xr = -sin(mu)*xe + cos(phi)*ye + sin(phi)*cos(mu)*ze;
yr = (sin(beta1)*sin(mu)+cos(alpha1)*cos(beta1)*sin(phi)*cos(mu))*ye -...
    (sin(alpha1)*sin(mu)+cos(alpha1)*cos(phi)*cos(mu))*ze;
zr = -cos(alpha1)*cos(beta1)*xe - sin(alpha1)*cos(beta1)*ye - sin(beta1)*ze;
%%�����Ӳ�ͨ���״�����ϵ�õ��ع�����ϵ,�Ӳ�����ECEF
coor_c = Rs*(sin(El)*cos(Az)*xr + sin(El)*sin(Az)*yr + cos(El)*zr);
%%������״ﻷ�������붯Ŀ���⼼���о�����3.3���������
%%�ȼ���Ϊ��Բ����ƫ����e=0������ƫ��f=0��
e = 0;
f = 0;
a = 6373e3; % ���򳤰��ᣬ�����õ�Բ����뾶m����Ϊ��ʱ��ΪΪԲ����
p = sqrt(coor_c(1)^2+coor_c(2)^2);
r = sqrt(p^2+coor_c(3)^2);
u = atan(coor_c(3)/p*((1-f)+e^2*a/r));
alpha2 = atan(coor_c(2)/coor_c(1)); %�Ӳ��ľ���
t1 = coor_c(3)*(1-f)+e^2*a*sin(u)^3;
t2 = (1-f)*(p-e^2*a*cos(u)^3);
beta2 = atan(t1/t2); %�Ӳ���γ��
h = p*cos(alpha2);
end

