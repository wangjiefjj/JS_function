function [ delta ] = fun_Delta( H )
%FUN_DELTA 此处显示有关此函数的摘要
%   此处显示详细说明
%% 式4.61
Ve = 0.4651; %km/s
Vp = fun_Vp(H)./1000; %km/s
Re = 6373; %flat earth radius，km
delta = Ve/Vp*(1+H/Re);

end

