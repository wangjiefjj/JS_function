function [ R_CC,alpha0 ] = fun_CCIter2(X,R,R_KA )
%FUN_CC 此处显示有关此函数的摘要
%   此处显示详细说明
%%迭代求解，训练样本估计的协方差和先验协方差的线性组合，利用凸优化得到组合系数。
%%X:训练样本
%R,样本估计的协方差
%R_KA:先验协方差
% R = fun_SCMN(X);
R_CC = R;
[N,K] = size(X);
alpha = 0;
alpha0 =0;
for k = 1:1
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
    for i = 1:K
        Rk = X(:,i) * X(:,i)';
        t1 = t1 + norm(Rk - R_CC,'fro')^2 + real(trace(R_KA*R_CC' - R_KA*Rk' + R_CC*Rk'-R_CC*R_CC'));
        t2 = t2 + norm(R_KA - Rk,'fro')^2;
    end
    alpha0 = max(min(1,(t1/t2)),0);   
    R_CC_0 = R_CC;
    R_CC = alpha0 * R_KA +(1-alpha0) * R_CC_0;
    if abs(alpha0- alpha)<1e-2
        break;
    end
end
end


