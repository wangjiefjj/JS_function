function [ FootLength ] = fun_FootLength( H,EL_beam )
%FUN_FOOTLENGTH 此处显示有关此函数的摘要
%   此处显示详细说明
% EL_beam: 俯仰向波束宽度。rad
% H:天基雷达距星下点的距离,km
%计算足印长度
Re = 6373e3; %flat earth radius，m
GrazeT = fun_GrazeTH(H,EL_beam,1);
GrazeH = fun_GrazeTH(H,EL_beam,2);
FootLength = Re.*(GrazeH-GrazeT-EL_beam);

end

