function [ R,rho] = fun_SFP( X , opt)
%FUN_SFP 此处显示有关此函数的摘要
%   此处显示详细说明
%%shrinkage fixed Point 方法
if nargin == 1
    opt = 1; 
end
[N,M] = size(X);%N:样本维数，M：样本个数 
R_NSCM = fun_NSCMN(X); %%基本的fixed Point估计结果
if opt == 1 %%<Robust Shrinkage Estimation of High-Dimensional Covariance Matrices>(17)式
    rho = (N^2 - 1/N*trace(R_NSCM*R_NSCM'))/((N^2-N*M-M) + (M+(M-1)/N)*trace(R_NSCM*R_NSCM'));
elseif opt == 2 %%<Regularized -Estimators of Scatter Matrix>(20)式
    %%或者《Optimal Design of the Adaptive Normalized Matched Filter Detector Using 
    %%Regularized Tyler Estimators》 （13）式
    t1 = N*trace(R_NSCM)-1;
    t2 = (N*trace(R_NSCM)-1) + M*(N+1)*( 1/N*trace(inv(R_NSCM*R_NSCM))-1 );
    rho = t1/t2;
end

R = eye(N,N);
for iter = 1:1
    Rt = zeros(N,N);
    R0 = R;
    %%fixed Point 方法<<Generalized Robust Shrinkage Estimator and Its
    % Application to STAP Detection Problem>> 式（3）
    for i = 1:M
        Rt = Rt + N/M*(X(:,i)*X(:,i)')/(X(:,i)'*inv(R0)*X(:,i)); 
    end
    %%%%%%
    Rt =(1-rho) * Rt + rho*eye(N,N);%
    R = N/trace(Rt)*Rt;
    if norm(R-R0)/norm(R)<1e-2
        break;
    end
end

end

