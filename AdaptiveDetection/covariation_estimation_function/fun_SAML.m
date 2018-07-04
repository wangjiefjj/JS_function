function [ R_SAML ] = fun_SAML( X,beta )
%%《Generalized Robust Shrinkage Estimator and Its
%Application to STAP Detection Problem》
%%<<Optimal Design of the Adaptive Normalized Matche
% Filter Detector Using Regularized Tyler Estimators>>
%   此处显示详细说明
%  收缩固定点方法，就是AML的对角加载
%X:训练样本
%beta: 对角加载量
%%Optimal Designd 的（13）式计算得来
if nargin == 1
    R_AML = fun_AML(X);
    [M,N] = size(X);
    t1 = (M*trace(R_AML)-1);
    t2 = M*trace(R_AML)-1+N*(M+1)*((1/M)*trace(R_AML^(2))-1);    
    beta = abs(t1/t2);
end

[M,N] = size(X);
I = eye(M,M);
R_SAML = eye(M,M);%R_x0;%eye(N,N);%以单位阵为迭代初值是，第二次迭代结果为NSCM结果
tao_child = 1;%%本次迭代值
tao_parent = 0;%%上次迭代值
count = 0;
while (abs(tao_child-tao_parent)>0.01)%
    tao_parent = tao_child;
    iR_SAML = inv(R_SAML);
    tao_child = diag(abs(X'*iR_SAML*X)/M);
    R_SAML_t = 0;
    for i = 1:N
        R_SAML_t = R_SAML_t+X(:,i)*X(:,i)'/N/tao_child(i);
    end
    R_SAML = ((1-beta)*R_SAML_t + beta * I);
    count = count + 1;
    if count > 6
        break;
    end
end
end

