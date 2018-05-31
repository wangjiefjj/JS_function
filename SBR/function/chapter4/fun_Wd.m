function [ wd ] = fun_Wd( H,R,Az, Tr, lambda)
%FUN_WD 此处显示有关此函数的摘要
%   此处显示详细说明
%% 无地球自转的多普勒频率
Re = 6373e3; %flat earth radius，m
Vp = fun_Vp(H*10^3);
% Vp = 7.61e3; %m/s
beta = 2*Vp*Tr/2/lambda;
t1 = Re*sin(R./Re).*cos(Az);
t2 = fun_R2Rs(H,R);
% t2 = sqrt(Re^2+(Re+H)^2-2*Re*(Re+H)*cos(R/Re));%%Rs
wd = beta *t1./t2;
end
