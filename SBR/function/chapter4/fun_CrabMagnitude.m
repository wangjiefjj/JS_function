function [ CrabM ] = fun_CrabMagnitude( alpha1, eta, H)
%FUN_CRABMAGNITUDE 此处显示有关此函数的摘要
%   此处显示详细说明
%% 偏航角计算
%% 参数说明
% alpha1: 雷达信息点的纬度
% eta：雷达轨道与赤道夹角
% H: 雷达距离星下点高度
%% 计算delta 4.61
Re = 6373e3; %flat earth radius，m
Ve = 0.4651e3; %赤道速度，m/s
Vp = fun_Vp(H); %m/s
delta = Ve./Vp.*(1+H./Re);
%% 计算偏航幅度 4.73
CrabM = sqrt(1+delta.^2.*cos(alpha1).^2-2.*delta.*cos(eta));
end

