function [ R_S ] = fun_SCMN( X,opt )
%����Э�������/N
%һ����һ�����뵥Ԫ
%X:ѵ������
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

