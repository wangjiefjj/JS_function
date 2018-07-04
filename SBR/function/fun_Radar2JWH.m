function [ alpha2, beta2, h] = fun_Radar2JWH( Az,El,Rs, alpha1,beta1,phi )
%FUN_RADAR2JWH 此处显示有关此函数的摘要
%   此处显示详细说明
% 知道波束照射的杂波块在SBR坐标系中的坐标，求杂波块在地球上的经纬高信息
% 根据《天基雷达空时自适应杂波抑制技术》P23、P24内容
%和《天基雷达环境仿真与动目标检测技术研究》（3.3）得到 
%% 输入 
% Az：波束方位角deg
% El：波束俯仰角deg
% Rs: 感兴趣点与SBR的斜距m
% alpha1: SBR的经度deg,东经为正，西经为负
% beta1： SBR的纬度deg ，北纬为正，南纬为负
% phi： SBR的轨道倾角deg
%% 开始
Az = Az/180*pi;
El = El/180*pi;
alpha1 = alpha1/180*pi;
beta1 = beta1/180*pi;
phi = phi/180*pi;
mu = asin(sin(beta1)/sin(phi));
%%假设地固坐标系的三个轴分别为
xe = [1,0,0];
ye = [0,1,0];
ze = [0,0,1];
%%则雷达的坐标系为
xr = -sin(mu)*xe + cos(phi)*ye + sin(phi)*cos(mu)*ze;
yr = (sin(beta1)*sin(mu)+cos(alpha1)*cos(beta1)*sin(phi)*cos(mu))*ye -...
    (sin(alpha1)*sin(mu)+cos(alpha1)*cos(phi)*cos(mu))*ze;
zr = -cos(alpha1)*cos(beta1)*xe - sin(alpha1)*cos(beta1)*ye - sin(beta1)*ze;
%%所以杂波通过雷达坐标系得到地固坐标系,杂波坐标ECEF
coor_c = Rs*(sin(El)*cos(Az)*xr + sin(El)*sin(Az)*yr + cos(El)*zr);
%%《天基雷达环境仿真与动目标检测技术研究》（3.3）迭代求解
%%先假设为正圆地球，偏心率e=0；椭球反偏率f=0；
e = 0;
f = 0;
a = 6373e3; % 椭球长半轴，这里用的圆地球半径m，因为暂时认为为圆地球
p = sqrt(coor_c(1)^2+coor_c(2)^2);
r = sqrt(p^2+coor_c(3)^2);
u = atan(coor_c(3)/p*((1-f)+e^2*a/r));
alpha2 = atan(coor_c(2)/coor_c(1)); %杂波的经度
t1 = coor_c(3)*(1-f)+e^2*a*sin(u)^3;
t2 = (1-f)*(p-e^2*a*cos(u)^3);
beta2 = atan(t1/t2); %杂波的纬度
h = p*cos(alpha2);
end

