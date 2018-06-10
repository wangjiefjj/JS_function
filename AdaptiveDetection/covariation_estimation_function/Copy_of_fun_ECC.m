function [ R_ECC,alpha ] = fun_ECC( X,R,R_KA,opt )
%FUN_ECC 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%%%%%%通过欧几里得几何均值得到的CC的权重因子估计
if nargin<4
    opt=1;
end
[N,K] = size(X);
t1 = 0;
t2 = 0;
if opt==1
    for i = 1:K
%         Rk = X(:,i) * X(:,i)';
        Rk = fun_Positive(X(:,i),4);
        t1 = t1 + norm(Rk - R,'fro')^2/K ;%+ real(trace(R_KA*R' - R_KA*Rk' + R*Rk'-R*R'))/K
        t2 = t2 + norm(R_KA - Rk,'fro')^2/K;
    end
else
    for i = 1:K
        Rk = X(:,i) * X(:,i)'/( X(:,i)'* X(:,i)/N);
        t1 = t1 + norm(Rk - R,'fro')^2/K + real(trace(R_KA*R' - R_KA*Rk' + R*Rk'-R*R'))/K;
        t2 = t2 + norm(R_KA - Rk,'fro')^2/K;
    end
    
end
alpha = max(min(1,(t1/t2)),0);   
R_ECC = alpha * R_KA +(1-alpha) * R;
end

