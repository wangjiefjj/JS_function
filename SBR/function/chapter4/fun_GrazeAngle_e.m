function [ graze ] = fun_GrazeAngle_e( H, R, Rs,alpha2)
%FUN_ 此处显示有关此函数的摘要
%   此处显示详细说明
%% 掠地角与地距的关系椭圆地球，公式（4B.12）(4B.34)
% H:天基雷达距星下点的距离,m
% R:地距m
% Rs:斜距m
% alpha2：目标纬度deg
%%
e = 0.08199;%地球离心率 
graze = fun_GrazeAngle(H,R,Rs); %圆地球掠射角deg
v = atan(e^2*sin(alpha2)/(2*(1-e^2*cos(alpha2)^2)))/pi*180; %修正量deg
graze = graze + v;
end

