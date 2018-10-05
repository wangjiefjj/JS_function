function [Trao] = fun_GRE_Rao(Train,x0,p)
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
x_ba = S^(-0.5)*x0;
p_ba = S^(-0.5)*p;
P_p = p_ba / (p_ba'*p_ba) * p_ba';
P_pp = eye(N)-P_p;
Trao = abs(x_ba'*P_p*x_ba)/abs((1+x_ba'*x_ba)*(1+x_ba'*P_pp*x_ba));
end

