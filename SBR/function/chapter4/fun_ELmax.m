function [ ELmax ] = fun_ELmax( H )
%FUN_ELMAX 此处显示有关此函数的摘要
%   此处显示详细说明
%根据H的最大的俯仰角
%%
Re = 6373; %flat earth radius，km
ELmax = pi/2 - (acos(1./(1+H/Re))); %公式（4.14）
ELmax = ELmax/pi*180;
end

