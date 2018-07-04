function [ F ] = fun_F( N_az, N_el, d_az, d_el, lambda, Az, El, Az0, El0, EL, Iaz,Iel)
%FUN_F 此处显示有关此函数的摘要
%   此处显示详细说明
%%方向图因子求解
%% 参数说明
% N_az: 方位向阵元数
% N_el: 俯仰向阵元数
% d_az: 方位向阵元间隔m
% d_el: 俯仰向阵元间隔m
% lambda: 波长
% Az: 要求的方向图方位角deg
% El: 要求的方向图俯仰角deg
% Az0: 主波束方位指向deg
% El0: 主波束俯仰指向deg
% EL:  天线法线和雷达坐标系Z轴所在平面夹角(天线倾斜角)
% Iaz: 方位向加权因子
% Iel: 俯仰向向加权因子
%% 方向图
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

