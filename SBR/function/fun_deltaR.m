function [ deltaR ] = fun_deltaR( H, R, Rs, Tr)
%FUN_DR 此处显示有关此函数的摘要
%   此处显示详细说明
% H: 卫星高度 km
% R: 地距 km
% Rs: 斜距 km
% dR: km
c = 3e8;  %光速m/s
%% 最大不模糊距离
Graze = fun_GrazeAngle(H,R,Rs)/180*pi;
deltaR =  c*Tr/2 *sec(Graze);
end

