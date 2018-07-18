function [ Ax ] = fun_Ax( N_az, N_el, d_az, d_el, lambda, Az, El, Az0, El0, Iaz,Iel )
%FUN_F 此处显示有关此函数的摘要
%   根据SBR书的(6.16)
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
% Iaz: 方位向加权因子
% Iel: 俯仰向向加权因子
%% 方向图
if nargin == 9 %%方位俯仰向都是均匀加权
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

