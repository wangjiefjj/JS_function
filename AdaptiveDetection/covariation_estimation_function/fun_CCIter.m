function [ R_CC,alpha0 ] = fun_CCIter(X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%������⣬ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
%R_KA:����Э����
% R = fun_SCMN(X);
[R_CC,alpha0_1] = fun_CC(X,R,R_KA);
for k = 1:10
    [R_CC,alpha0] = fun_CC(X,R_CC,R_KA);
    if abs(alpha0_1 - alpha0) <1e-2
        break;
    end
    alpha0_1 = alpha0;
end
end


