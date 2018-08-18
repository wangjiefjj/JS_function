function [ R_CC,alpha0 ] = fun_CCIter3(X,R,R_KA )
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
Rkk = zeros(N,N);
for k = 1:1
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
    for i = 1:K
        Rk = alpha0 * R_KA +(1-alpha0) *(X(:,i) * X(:,i)');
        Rkk = Rkk+ Rk/K;
        t1 = t1 + norm(Rk - R_CC,'fro')^2 ;
        t2 = t2 + norm(R_KA - Rk,'fro')^2;
    end
    alpha0 = (t1/t2);
%     alpha0 = alpha0^3;
    R_CC = alpha0 * R_KA +(1-alpha0) * Rkk;
    Rkk = zeros(N,N);
%     if abs(alpha0- alpha)<1e-3
%         break;
%     end
end
end


