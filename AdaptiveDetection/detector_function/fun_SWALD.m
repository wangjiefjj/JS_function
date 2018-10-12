function [Tswald] = fun_SWALD(Train,x0,H)
%FUN_SWALD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%�ӿռ�Wald
S = Train*Train';
x = S^(-0.5)*x0;
H = S^(-0.5)*H;
P_H = H/(H'*H)*H';
Tswald = abs(x'*P_H*x);
end

