function [All] = fun_ComputeAzEl( alpha1,beta1,h1,eta,alpha2,beta2,h2 )
%FUN_COMPUTEAZEL 此处显示有关此函数的摘要
%   此处显示详细说明
%  根据目标经纬高和雷达经纬高轨道倾角来计算波束的指向
%alpha1: 雷达经度
%beta1: 雷达纬度
%h1: 雷达高（H）
%eta: 雷达轨道倾角
Re = 6373e3;
[X1,Y1,Z1] = fun_JWH2XYZ(alpha1,beta1, h1+Re); %卫星经纬高转XYZ
[X2,Y2,Z2] = fun_JWH2XYZ(alpha2,beta2, h2+Re);   %目标点经纬高转XYZ
All.Rs = sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);
All.R = fun_Rs2R(h1,All.Rs);
%计算波束指向的俯仰角
All.graze = fun_GrazeAngle(h1,All.R,All.Rs);
All.El = fun_ELAngle(h1,All.R);
% 卫星雷达坐标系坐标
All.mu = asin(sin(beta1/180*pi)/sin(eta/180*pi));
[All.Xr1,All.Yr1,All.Zr1] = fun_XYZ2Radar(alpha1,beta1,eta,Re+h1, X1,Y1,Z1);
% Xr1 = -sin(mu)*X1+cos(eta/180*pi)*cos(mu)*Y1+sin(eta/180*pi)*cos(mu)*(Z1-Re-H);
All.Az = abs(acos(All.Xr1/All.Rs/sin(All.El/180*pi)))/pi*180;
end

