function [ R ] = fun_RPowerEMean( X,alpha,opt)
%FUN_RPOWEREMEAN 此处显示有关此函数的摘要
%   此处显示详细说明
%%Power-E的几何均值协方差
if nargin ==1
    alpha = 2;
end
[N,L] = size(X);
R_alpha = zeros(N,N);
for i = 1:L
     Ri = fun_Positive(X(:,i),opt);
     [UA, LA] = svd(Ri);
     Ri = UA * diag((diag(LA).^alpha)) * UA';
%      Ri = Ri^alpha;
     R_alpha = R_alpha + Ri/L;
end
alpha = 1/alpha;
[UA, LA] = svd(R_alpha);
R = UA * diag((diag(LA).^alpha)) * UA';
% R = R_alpha^alpha;


end

