function [Tsglrt] = fun_SGLRT(Train,x0,H)
%FUN_SGLRT �˴���ʾ�йش˺�����ժҪ
%�ӿռ�GLRT
S = Train*Train';
x = S^(-0.5)*x0;
H = S^(-0.5)*H;
P_H = H/(H'*H)*H';
Tsglrt = abs(x'*P_H*x)/abs(1+x'*x-x'*P_H*x);
end

