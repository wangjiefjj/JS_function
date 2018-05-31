function [ GA,R,Rs ] = fun_GrazeAngle( H, R, Rs)
%FUN_ 此处显示有关此函数的摘要
%   此处显示详细说明
%% 掠地角与地距的关系，公式（4.6）
% H:天基雷达距星下点的距离,km
%%
Re = 6373; %flat earth radius，km
if nargin==1
    [R,Rs] = fun_RsR(H);
end
GA = abs(acos((Re+H)./Rs .* sin(R/Re)));
GA = GA/pi*180;
end

