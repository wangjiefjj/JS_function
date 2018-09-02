function [ R ] = fun_FP( X )
%FUN_FP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%fixed Point ����<<Generalized Robust Shrinkage Estimator and Its
% Application to STAP Detection Problem>> ʽ��3��
% X��ѵ������
[m,N] = size(X);%m:ѵ����������ά����Nѵ��������������
R = eye(m,m);
for iter = 1:10
    Rt = zeros(m,m);
    R0 = R;
    for i = 1:N
        %%��һ��
        %%<Robust Shrinkage Estimation of High-DimensionalCovariance Matrices>
        %%��4��ʽ
%         X(:,i) = X(:,i)/sqrt(norm(X(:,i),'fro')^2/m);
%         X(:,i) = X(:,i)/norm(X(:,i));
        Rt = Rt + m/N*(X(:,i)*X(:,i)')/(X(:,i)'*inv(R0)*X(:,i)); 
    end
    R = m/trace(Rt)*Rt;
    if norm(R-R0)/norm(R)<1e-2
        break;
    end
end
end

