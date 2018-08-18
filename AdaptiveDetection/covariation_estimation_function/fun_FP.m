function [ R ] = fun_FP( X )
%FUN_FP 此处显示有关此函数的摘要
%   此处显示详细说明
%%fixed Point 方法<<Generalized Robust Shrinkage Estimator and Its
% Application to STAP Detection Problem>> 式（3）
% X：训练数据
[m,N] = size(X);%m:训练数据向量维数，N训练数据向量个数
R0 = eye(m,m);
R = eye(m,m);
for iter = 1:10
    Rt = zeros(m,m);
    R0 = R;
    for i = 1:N
        Rt = Rt + m/M*(X(:,i)*X(:,i)')/(X(:,i)'*R0*X(:,i)); 
    end
    R = m/trace(Rt)*Rt;
    if norm(R-R0)/norm(R)<1e-2
        break;
    end
end
end

