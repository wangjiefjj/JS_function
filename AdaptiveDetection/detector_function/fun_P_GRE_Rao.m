function [Tprao] = fun_P_GRE_Rao(Train,x0,H)
%FUN_PGRE 此处显示有关此函数的摘要
%   此处显示详细说明
%%GRE关系下的Persymmetric  Rao 检测器 
[N,L] = size(Train);
%%%产生置换矩阵，反对角线都为1其余为0
J = zeros(N,N);
for i = 1:N
    J(i,N-i+1) = 1;
end
S = Train*Train';
Sp = 0.5*(S + J*conj(S)*J);
Y = [0.5*(x0+J*conj(x0)),0.5*(x0-J*conj(x0))];
S0 = Sp + Y*Y';
x_ba = S0^(-0.5)*x0;
p_ba = S0^(-0.5)*H;
P_p = p_ba / (p_ba'*p_ba) * p_ba';
Tprao = abs(x_ba'*P_p*x_ba);
end

