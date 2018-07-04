function [ F ] = fun_F( N_az, N_el, d_az, d_el, lambda, Az, El, Az0, El0, EL, Iaz,Iel)
%FUN_F �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%����ͼ�������
%% ����˵��
% N_az: ��λ����Ԫ��
% N_el: ��������Ԫ��
% d_az: ��λ����Ԫ���m
% d_el: ��������Ԫ���m
% lambda: ����
% Az: Ҫ��ķ���ͼ��λ��deg
% El: Ҫ��ķ���ͼ������deg
% Az0: ��������λָ��deg
% El0: ����������ָ��deg
% EL:  ���߷��ߺ��״�����ϵZ������ƽ��н�(������б��)
% Iaz: ��λ���Ȩ����
% Iel: ���������Ȩ����
%% ����ͼ
if nargin == 9
    Iaz = ones(1,N_az);
    Iel = ones(1,N_el);
end
SinAz = sin(Az/180*pi);
CosAz = cos(Az/180*pi);
SinEl = sin(El/180*pi);
CosEl = cos(El/180*pi);
SinEl0 = sin(El0/180*pi);
CosEl0 = cos(El0/180*pi);
SinAz0 = sin(Az0/180*pi);
CosAz0 = cos(Az0/180*pi);
SinEL = sin(EL/180*pi);
CosEL = cos(EL/180*pi);
kx = 2*pi*d_az/lambda;
ky = 2*pi*d_el/lambda;
tx = SinEl*CosAz - SinEl0*CosAz0;
ty = CosEL*(SinEl*SinAz-SinEl0*SinAz0) - SinEL*(CosEl-CosEl0);
t1 = Iaz.*exp(1j*(0:N_az-1)*kx*tx);
t2 = Iel.*exp(1j*(0:N_el-1)*ky*ty);
F = sum(kron(t1,t2));
end

