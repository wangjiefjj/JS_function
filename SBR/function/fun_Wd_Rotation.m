function [ wd_r ] = fun_Wd_Rotation( H,R,Az, Tr, lambda,alpha1,eta,~)
%FUN_WD 此处显示有关此函数的摘要
%   此处显示详细说明
%Az: 方位角rad
%% 有地球自转的多普勒频率
if nargin ==7
    opt=1;%目标在卫星轨道东边
else
    opt=-1;%目标在卫星轨道西边边
end
%% 
Re = 6373e3; %flat earth radius，m
Vp = fun_Vp(H*10^3);
% Vp = 7.61e3; %m/s
beta = 2*Vp*Tr/2/lambda;
crabA = abs(fun_CrabAngle(alpha1,eta,H));
crabM = abs(fun_CrabMagnitude(alpha1,eta,H));
t1 = Re*sin(R./Re).*cos(Az+opt*crabA);
t2 = fun_R2Rs(H,R);
% t2 = sqrt(Re^2+(Re+H)^2-2*Re*(Re+H)*cos(R/Re));%%Rs
wd_r = beta .*crabM.*t1./t2;


end

