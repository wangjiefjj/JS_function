function [ R_CC,alpha0 ] = fun_CCIter3(X,R,R_KA )
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
Rkk = zeros(N,N);
for k = 1:1
    alpha = alpha0;
    t1 = 0;
    t2 = 0;
    for i = 1:K
        Rk = alpha0 * R_KA +(1-alpha0) *(X(:,i) * X(:,i)');
        Rkk = Rkk+ Rk/K;
        t1 = t1 + norm(Rk - R_CC,'fro')^2 ;
        t2 = t2 + norm(R_KA - Rk,'fro')^2;
    end
    alpha0 = (t1/t2);
%     alpha0 = alpha0^3;
    R_CC = alpha0 * R_KA +(1-alpha0) * Rkk;
    Rkk = zeros(N,N);
%     if abs(alpha0- alpha)<1e-3
%         break;
%     end
end
end


