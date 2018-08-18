function [ R,alpha ] = fun_LogCC_new( X,R_KA )
%FUN_LOGCC_NEW 此处显示有关此函数的摘要
%   此处显示详细说明
%logE推导的色加载
[N,L] = size(X);
Ri = zeros(N,N,L);
for i = 1:L
    Ri(:,:,i) = fun_Positive(X(:,i),4);
%     logm_R = logm_R + fun_Logm(Ri(:,:,i))/L;
end
t1 = 0;
t2 = 0;
logm_R = logm(fun_RLogEMean(X,1));
for i = 1:L
    t11 = logm(Ri(:,:,i)) - logm_R;
    t1 = t1 + norm(t11,'fro')^2;
    t22 = logm(R_KA) - logm(Ri(:,:,i));
    t2 = t2 + norm(t22,'fro')^2;
end
alpha = min(1,t1/t2);
R = alpha * logm(R_KA) + (1-alpha)*logm_R;
R = expm(R);