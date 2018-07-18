function [ Xr,Yr,Zr ] = fun_XYZ2Radar( alpha1,beta1,eta,h,Xe,Ye,Ze)
%FUN_XYZ2RADAR 此处显示有关此函数的摘要
%   此处显示详细说明
%%地心坐标系转雷达坐标系
% alpha1: SBR的经度deg,东经为正，西经为负
% beta1： SBR的纬度deg ，北纬为正，南纬为负
% eta： SBR的轨道倾角deg
% h 雷达高度（Re+H）
% Xe:雷达在地固坐标系的坐标
alpha1 = alpha1/180*pi;
beta1 = beta1/180*pi;
eta = eta/180*pi;
mu = asin(sin(beta1)/sin(eta));
Ar2e = [-sin(mu),         0,       -cos(alpha1)*cos(beta1);
        cos(eta)*cos(mu), sin(beta1)*sin(mu)+cos(alpha1)*cos(beta1)*sin(eta)*cos(mu),  -sin(alpha1)*cos(beta1);
        sin(eta)*cos(mu), -(sin(alpha1)*sin(mu)+cos(alpha1)*cos(eta)*cos(mu)), - sin(beta1)];
Ae2r = Ar2e.';
Z = Ae2r * [Xe;Ye;Ze-h];
Xr = Z(1);
Yr = Z(2);
Zr = Z(3);
end
