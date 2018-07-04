function [ sigma0 ] = fun_Morchin( graze, fo, opt, ss)
%FUN_MOR 此处显示有关此函数的摘要
%   此处显示详细说明
%Morchin模型，模拟地杂波
%<<地杂波背景中机载预警雷达作用距离分析>>
% graze: 掠地角
% fo: 载频
% opt: 地形选择
% ss: 海情等级 0~9
c = 3e8;
lambda = c/fo;
mu = sqrt(fo/10^9)/4.7;
if nargin == 3
    ss = 1;
end
if opt == 1 %%沙漠
    A = 0.00126;
    B = pi/2;
    beta0 = 0.14;
    phic = asin(lambda/4/pi/9.3/(beta0^2.2));
    if graze<phic
       sigmac0 = graze/phic;
    else
       sigmac0 = 1;
    end
    
elseif opt == 2 %%农田
    A = 0.004;
    B = pi/2;
    beta0 = 0.2;
    sigmac0 = 1;
elseif opt == 3 %%丘陵
    A = 0.0126;
    B = pi/2;
    beta0 = 0.4;
    sigmac0 = 1;
elseif opt == 4 %%高山
    A = 0.04;
    B = 1.24;
    beta0 = 0.5;
    sigmac0 = 1;
elseif opt == 5 %%海洋
    %<地(海)杂波反射率模型研究,彭世蕤,汤子跃>
    A = 4e-7*10^(0.6*(ss+1));
    B = pi/2;
    beta0 = 2.44*((ss+1)^(1.08))/57.29;
    he = 0.025+0.046*ss^(1.72);
    phic = asin(lambda/4/pi/he);
    k = 1.9;
    if graze<phic
        sigmac0 = (graze/phic)^k;
    else
        sigmac0 = 1;
    end
    mu = 1;
end 
sigma0 = A*sigmac0*sin(graze)/lambda + ...
    mu*cot(beta0)^2*exp(-tan(B-graze)^2/tan(beta0)^2);
end

