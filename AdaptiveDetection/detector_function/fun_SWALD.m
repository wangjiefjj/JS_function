function [Tswald] = fun_SWALD(Train,x0,H)
%FUN_SWALD 此处显示有关此函数的摘要
%   此处显示详细说明
%子空间Wald
S = Train*Train';
x = S^(-0.5)*x0;
H = S^(-0.5)*H;
P_H = H/(H'*H)*H';
Tswald = abs(x'*P_H*x);
end

