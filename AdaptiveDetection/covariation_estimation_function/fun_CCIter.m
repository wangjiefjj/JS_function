function [ R_CC,alpha0 ] = fun_CCIter(X,R,R_KA )
%FUN_CC 此处显示有关此函数的摘要
%   此处显示详细说明
%%迭代求解，训练样本估计的协方差和先验协方差的线性组合，利用凸优化得到组合系数。
%%X:训练样本
%R,样本估计的协方差
%R_KA:先验协方差
% R = fun_SCMN(X);
[R_CC,alpha0_1] = fun_CC(X,R,R_KA);
for k = 1:10
    [R_CC,alpha0] = fun_CC(X,R_CC,R_KA);
    if abs(alpha0_1 - alpha0) <1e-2
        break;
    end
    alpha0_1 = alpha0;
end
end


