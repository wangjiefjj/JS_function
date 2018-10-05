function [Tprao] = fun_P_GRE_Rao(Train,x0,H)
%FUN_PGRE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%GRE��ϵ�µ�Persymmetric  Rao ����� 
[N,L] = size(Train);
%%%�����û����󣬷��Խ��߶�Ϊ1����Ϊ0
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

