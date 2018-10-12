function [Tsrao] = fun_SRAO(Train,x0,H)
%FUN_SGLRT �˴���ʾ�йش˺�����ժҪ
%�ӿռ�RAO
[N,~]=size(Train);
S = Train*Train';
x = S^(-0.5)*x0;
H = S^(-0.5)*H;
P_H = H/(H'*H)*H';
P_Hp = eye(N) - P_H;
Tsrao = abs(x'*P_H*x)/abs((1+x'*x)*(1+x'*P_Hp*x));
end

