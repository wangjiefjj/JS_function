function [Tpwald] = fun_P_GRE_Wald(Train,x0,H)
%FUN_P_GRE_WALD 此处显示有关此函数的摘要
%   此处显示详细说明
%%GRE关系下的Persymmetric  Wald 检测器 
[N,L] = size(Train);
%%%产生置换矩阵，反对角线都为1其余为0
J = zeros(N,N);
for i = 1:N
    J(i,N-i+1) = 1;
end
S = Train*Train';
Sp = 0.5*(S + J*conj(S)*J);
% Y = [0.5*(x0+J*conj(x0)),0.5*(x0-J*conj(x0))];
% Theta = (H'/Sp*H)\(H'/Sp*x0);
% Thetap = (Y'/Sp*H)/(H'/Sp*H);
% % % Theta = Thetap(:,1)+1j*Thetap(:,2);
% S1 = Sp + (Y-H*Thetap)*(Y-H*Thetap)';
% vt = H*Theta;
% Tpwald = abs(vt'/S1*vt);

x = Sp^(-0.5)*x0;
H = Sp^(-0.5)*H;
S = Sp^(-0.5)*S*Sp^(-0.5);
P_H = H/(H'*H)*H';
iS = inv(S);
Tpwald = abs(x'*P_H*x);
end

