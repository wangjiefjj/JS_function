function [ R_SAML ] = fun_SAML( X,beta )
%%��Generalized Robust Shrinkage Estimator and Its
%Application to STAP Detection Problem��
%%<<Optimal Design of the Adaptive Normalized Matche
% Filter Detector Using Regularized Tyler Estimators>>
%   �˴���ʾ��ϸ˵��
%  �����̶��㷽��������AML�ĶԽǼ���
%X:ѵ������
%beta: �ԽǼ�����
%%Optimal Designd �ģ�13��ʽ�������
if nargin == 1
    R_AML = fun_AML(X);
    [M,N] = size(X);
    t1 = (M*trace(R_AML)-1);
    t2 = M*trace(R_AML)-1+N*(M+1)*((1/M)*trace(R_AML^(2))-1);    
    beta = abs(t1/t2);
end

[M,N] = size(X);
I = eye(M,M);
R_SAML = eye(M,M);%R_x0;%eye(N,N);%�Ե�λ��Ϊ������ֵ�ǣ��ڶ��ε������ΪNSCM���
tao_child = 1;%%���ε���ֵ
tao_parent = 0;%%�ϴε���ֵ
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

