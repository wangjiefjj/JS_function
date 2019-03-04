function [s] = fun_RectSteer(N_az, N_el, d_az, d_el, lambda,Az0, El0,opt)
%FUN_RECTSTEER 此处显示有关此函数的摘要
%   此处显示详细说明
%% 矩形阵列的导向矢量形成
%% 参数说明
% N_az: 方位向阵元数
% N_el: 俯仰向阵元数
% d_az: 方位向阵元间隔m
% d_el: 俯仰向阵元间隔m
% lambda: 波长
% Az0: 导向矢量方位角指向deg
% El0: 导向矢量俯仰指向deg
% opt: 选项，1：s列向量，2：s为矩阵，默认opt=1
%% 导向矢量形成
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

