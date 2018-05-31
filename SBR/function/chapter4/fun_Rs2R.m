function [ R ] = fun_Rs2R( H, Rs )
%FUN_RS2R 此处显示有关此函数的摘要
%   此处显示详细说明
%% 指定高度和斜距求对应的地距
Re = 6373; %flat earth radius，km
t1 = Re^2+(Re+H)^2 - Rs.^2;
t2 = 2*Re*(Re+H);
R = acos(t1./t2).*Re;

end

