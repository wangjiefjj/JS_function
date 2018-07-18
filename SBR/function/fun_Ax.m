function [ Ax ] = fun_Ax( N_az, N_el, d_az, d_el, lambda, Az, El, Az0, El0, Iaz,Iel )
%FUN_F �˴���ʾ�йش˺�����ժҪ
%   ����SBR���(6.16)
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
% Iaz: ��λ���Ȩ����
% Iel: ���������Ȩ����
%% ����ͼ
if nargin == 9 %%��λ�������Ǿ��ȼ�Ȩ
    Iaz = ones(1,N_az);
    Iel = ones(1,N_el);
end
CosAz = cos(Az/180*pi);
CosEl = cos(El/180*pi);
CosEl0 = cos(El0/180*pi);
CosAz0 = cos(Az0/180*pi);
kx = 2*pi*d_az/lambda;
ky = 2*pi*d_el/lambda;
tx = CosAz - CosAz0;
ty = CosEl - CosEl0;
t1 = Iaz.*exp(1j*(0:N_az-1)*kx*tx);
t2 = Iel.*exp(1j*(0:N_el-1)*ky*ty);
Ax = sum(kron(t1,t2));

end

