function [ Rs ] = fun_R2Rs( H,R )
%FUN_RS2R 此处显示有关此函数的摘要
%   此处显示详细说明
%% 指定高度和地距求对应的斜距
Re = 6373; %flat earth radius，km
Rs = sqrt(Re^2 + (Re+H)^2 - 2*Re*(Re+H).*cos(R/Re));
end

