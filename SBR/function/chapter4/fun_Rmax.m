function [ Rmax ] = fun_Rmax( H )
%FUN_RMAX 此处显示有关此函数的摘要
%   此处显示详细说明
%高度和最大的地距
Re = 6373; %flat earth radius，km
Rmax = Re*acos(1./(1+H./Re));%公式（4.12）
end

