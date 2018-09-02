function [ C ] = fun_corrcoef( X )
%FUN_CORRCOEF 此处显示有关此函数的摘要
%   此处显示详细说明
%%相关系数矩阵，计算方法为华小强大论文，103页
N = length(X);
c = zeros(N,1);
for i = 0:N-1 % 行
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

