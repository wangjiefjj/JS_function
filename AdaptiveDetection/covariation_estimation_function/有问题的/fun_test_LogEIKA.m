function [ R_CC,alpha0 ] = fun_test_LogEIKA( X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%����log-E�����ĵ���CCЭ�������
%   �˴���ʾ��ϸ˵��
%%ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
% R = abs(fun_SCMN(X));
%R_KA:����Э����
[R_CC,alpha0_1] = fun_test_LogEKA(X,R,R_KA);
for k = 1:100
    [R_CC,alpha0] = fun_test_LogEKA(X,R_CC,R_KA);
    if abs(alpha0_1 - alpha0) <1e-2
        break;
    end
    alpha0_1 = alpha0;
end
end


