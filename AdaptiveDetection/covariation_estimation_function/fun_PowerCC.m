function [ R,alpha ] = fun_PowerCC( X,R_KA,beta,opt )
%FUN_POWERCC 此处显示有关此函数的摘要
%   此处显示详细说明
%logE推导的色加载
%btea powerE 的幂次
%opt HPD方法
if nargin == 2
    beta = 2;
end
[N,L] = size(X);
R_beta = fun_RPowerEMean( X,beta,opt);
Ri = zeros(N,N,L);
for i = 1:L
    Ri(:,:,i) = fun_Positive(X(:,i),opt);
%     logm_R = logm_R + fun_Logm(Ri(:,:,i))/L;
end
t1 = 0;
t2 = 0;
for i = 1:L
    t11 = Ri(:,:,i)^beta - R_beta^beta;
    t1 = t1 + norm(t11,'fro')^2;
    t22 = R_KA^beta - Ri(:,:,i)^beta;
    t2 = t2 + norm(t22,'fro')^2;
end
alpha = min(1,t1/t2);
R = alpha * R_KA + (1-alpha)*R_beta;
end

