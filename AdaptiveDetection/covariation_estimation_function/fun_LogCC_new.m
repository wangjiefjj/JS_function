function [ R,alpha ] = fun_LogCC_new( X,R_KA,opt )
%FUN_LOGCC_NEW 此处显示有关此函数的摘要
%   此处显示详细说明
%logE推导的色加载
[N,L] = size(X);
Ri = zeros(N,N,L);
for i = 1:L
    Ri(:,:,i) = fun_Positive(X(:,i),opt);
%     logm_R = logm_R + fun_Logm(Ri(:,:,i))/L;
end
t1 = 0;
t2 = 0;
if opt == 4 
    logm_R = logm(fun_RLogEMean(X,opt));
else
    logm_R = fun_Logm(fun_RLogEMean(X,opt)); 
end

if opt == 4
    for i = 1:L
        t11 = logm(Ri(:,:,i)) - logm_R;
        t1 = t1 + norm(t11,'fro')^2;
        t22 = logm(R_KA) - logm(Ri(:,:,i));
        t2 = t2 + norm(t22,'fro')^2;
    end
else
   for i = 1:L
        t11 = fun_Logm(Ri(:,:,i)) - logm_R;
        t1 = t1 + norm(t11,'fro')^2;
        t22 = fun_Logm(R_KA) - logm(Ri(:,:,i));
        t2 = t2 + norm(t22,'fro')^2;
    end 
end

alpha = min(1,t1/t2);
R = alpha * logm(R_KA) + (1-alpha)*logm_R;
if opt == 4
    R = expm(R);
else
    R = fun_Expm(R);
end
