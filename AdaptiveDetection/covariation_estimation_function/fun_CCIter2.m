function [ R_CC,alpha0 ] = fun_CCIter2(X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%������⣬ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
%R_KA:����Э����
% R = fun_SCMN(X);
R_CC = R;
[N,K] = size(X);
alpha = 0;
alpha0 =0;
for k = 1:1
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
    for i = 1:K
        Rk = X(:,i) * X(:,i)';
        t1 = t1 + norm(Rk - R_CC,'fro')^2 + real(trace(R_KA*R_CC' - R_KA*Rk' + R_CC*Rk'-R_CC*R_CC'));
        t2 = t2 + norm(R_KA - Rk,'fro')^2;
    end
    alpha0 = max(min(1,(t1/t2)),0);   
    R_CC_0 = R_CC;
    R_CC = alpha0 * R_KA +(1-alpha0) * R_CC_0;
    if abs(alpha0- alpha)<1e-2
        break;
    end
end
end


