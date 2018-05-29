function [ R,Rs ] = fun_RsR( H )
%FUN_RSR 此处显示有关此函数的摘要
%   此处显示详细说明：
%% 协矩Rs和地距R的关系，公式（4.3）
% H:天基雷达距星下点的距离,km
%%
Re = 6373; %flat earth radius，km
R = 0:fun_Rmax(H);
Rs = sqrt(Re^2 + (Re+H)^2 - 2*Re*(Re+H)*cos(R/Re));


end

