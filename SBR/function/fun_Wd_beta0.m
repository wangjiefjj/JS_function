function [ wd ] = fun_Wd_beta0( H,R,Az,beta)
%FUN_WD 此处显示有关此函数的摘要
%   此处显示详细说明
%% 无地球自转的多普勒频率
Re = 6373e3; %flat earth radius，m
t1 = Re*sin(R./Re).*cos(Az);
t2 = fun_R2Rs(H,R);
wd = beta *t1./t2;
end
