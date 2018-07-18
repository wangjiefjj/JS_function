function [EL] = fun_ELAngle_graze(H,graze)
%FUN_ELANGLE 此处显示有关此函数的摘要
%   此处显示详细说明
%% 俯仰角与地距的关系，公式（4.8）
%% H:高度m，graze：掠地角deg,EL:俯仰角deg
Re = 6373e3; %flat earth radius，m
EL = asin(1/(1+H/Re)*cos(graze/180*pi))/pi*180;

end

