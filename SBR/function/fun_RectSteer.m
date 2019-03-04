function [s] = fun_RectSteer(N_az, N_el, d_az, d_el, lambda,Az0, El0,opt)
%FUN_RECTSTEER �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �������еĵ���ʸ���γ�
%% ����˵��
% N_az: ��λ����Ԫ��
% N_el: ��������Ԫ��
% d_az: ��λ����Ԫ���m
% d_el: ��������Ԫ���m
% lambda: ����
% Az0: ����ʸ����λ��ָ��deg
% El0: ����ʸ������ָ��deg
% opt: ѡ�1��s��������2��sΪ����Ĭ��opt=1
%% ����ʸ���γ�
if nargin == 7 
    opt=1;
end
d_az = 2*pi*d_az/lambda;
d_el = 2*pi*d_el/lambda;
s_az = exp(1j*(0:N_az-1)*d_az*sin(deg2rad(El0))*cos(deg2rad(Az0)));
s_el = exp(1j*(0:N_el-1)*d_el*sin(deg2rad(El0))).';

s = s_el*s_az;
if opt==1
    s = reshape(s,N_az*N_el,1);
end

end

