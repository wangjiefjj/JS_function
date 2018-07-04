function [ Vp ] = fun_Vp( H )
%FUN_VP 此处显示有关此函数的摘要
%   此处显示详细说明
%H：高度单位为m
%% 根据卫星高度计算卫星速度
G = 6.673e-11; %%宇宙常数
Me = 5.965e24; %%地球质量kg
Re = 6373*10^3; %flat earth radius，m
Vp = sqrt(G*Me./(Re+H));%m/s
end

