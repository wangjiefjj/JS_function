function [ C ] = fun_corrcoef( X )
%FUN_CORRCOEF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%���ϵ�����󣬼��㷽��Ϊ��Сǿ�����ģ�103ҳ
N = length(X);
c = zeros(N,1);
for i = 0:N-1 % ��
    c(i+1) = sum(X(1:N-i) .* conj(X(i+1:N))) / (N-i);%(N-i)
end
% c = c.';
% c = conj(c); 
C = toeplitz(c);
for i = 1:N
    for j = i+1:N
        C(i,j) = conj(C(i,j));
    end
end
end

