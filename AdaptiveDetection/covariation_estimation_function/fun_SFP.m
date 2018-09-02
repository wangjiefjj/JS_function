function [ R,rho] = fun_SFP( X , opt)
%FUN_SFP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%shrinkage fixed Point ����
if nargin == 1
    opt = 1; 
end
[N,M] = size(X);%N:����ά����M���������� 
R_NSCM = fun_NSCMN(X); %%������fixed Point���ƽ��
if opt == 1 %%<Robust Shrinkage Estimation of High-Dimensional Covariance Matrices>(17)ʽ
    rho = (N^2 - 1/N*trace(R_NSCM*R_NSCM'))/((N^2-N*M-M) + (M+(M-1)/N)*trace(R_NSCM*R_NSCM'));
elseif opt == 2 %%<Regularized -Estimators of Scatter Matrix>(20)ʽ
    %%���ߡ�Optimal Design of the Adaptive Normalized Matched Filter Detector Using 
    %%Regularized Tyler Estimators�� ��13��ʽ
    t1 = N*trace(R_NSCM)-1;
    t2 = (N*trace(R_NSCM)-1) + M*(N+1)*( 1/N*trace(inv(R_NSCM*R_NSCM))-1 );
    rho = t1/t2;
end

R = eye(N,N);
for iter = 1:1
    Rt = zeros(N,N);
    R0 = R;
    %%fixed Point ����<<Generalized Robust Shrinkage Estimator and Its
    % Application to STAP Detection Problem>> ʽ��3��
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

