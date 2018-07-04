function [ Rsmax ] = fun_Rsmax( H )
%FUN_RSMAX 此处显示有关此函数的摘要
%   此处显示详细说明
%高度和最大的斜距
Re = 6373e3; %flat earth radius，m
Rsmax = sqrt(H.*(2*Re+H));%公式（4.13）

end

