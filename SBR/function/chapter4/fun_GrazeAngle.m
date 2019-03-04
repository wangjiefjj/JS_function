function [ GA,R,Rs ] = fun_GrazeAngle( H, R, Rs)
%FUN_ 此处显示有关此函数的摘要
%   此处显示详细说明
%% 掠地角与地距的关系，公式（4.6）
% H:天基雷达距星下点的距离,m
%%
Re = 6373e3; %flat earth radius，m
if nargin==1
    [R,Rs] = fun_RsR(H);
end
GA = (acos((Re+H)./Rs .* sin(R/Re)));
GA = rad2deg(GA);
end

