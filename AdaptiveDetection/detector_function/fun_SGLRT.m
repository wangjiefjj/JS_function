function [Tsglrt] = fun_SGLRT(Train,x0,H)
%FUN_SGLRT 此处显示有关此函数的摘要
%子空间GLRT
S = Train*Train';
x = S^(-0.5)*x0;
H = S^(-0.5)*H;
P_H = H/(H'*H)*H';
Tsglrt = (x'*P_H*x)/(1+x'*x-x'*P_H*x);
end

