function [ R_S ] = fun_SCMN( X,opt )
%采样协方差矩阵/N
%一列是一个距离单元
%X:训练样本
if nargin==1
    opt = 1;
end

[~,N] = size(X);
if opt==1 
    R_S = (X*X'/N);
else
    R_S = (X*X'/N);
end
end

