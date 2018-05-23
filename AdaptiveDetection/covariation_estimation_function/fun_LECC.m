function [ R_LECC,alpha ] = fun_LECC( X,R,R_KA,opt )
%FUN_ECC 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%%%%%%通过Log欧几里得几何均值得到的CC的权重因子估计
if nargin<4
    opt=1;
end
[N,K] = size(X);
t1 = 0;
t2 = 0;
R_KA = logm(R_KA);
R = logm(R);
if opt==1
    for i = 1:K
        Rk = logm(fun_Positive(X(:,i),4));
        t1 = t1 + norm(Rk - R,'fro')^2 + real(trace(R_KA*R' - R_KA*Rk' + R*Rk'-R*R'));
        t2 = t2 + norm(R_KA - Rk,'fro')^2/K;
    end
elseif opt==2
    for i = 1:K
        Rk = logm(fun_Positive(X(:,i),5));
        t1 = t1 + norm(Rk - R,'fro')^2 + real(trace(R_KA*R' - R_KA*Rk' + R*Rk'-R*R'));
        t2 = t2 + norm(R_KA - Rk,'fro')^2;
    end
    
end
alpha = max(min(1,(t1/t2)),0);   
R_LECC = alpha * R_KA +(1-alpha) * R;
R_LECC = expm(R_LECC);
end

