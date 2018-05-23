function [ R_CC,alpha0 ] = fun_test_LogEKA( X,R,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%����log-E������CCЭ�������
%   �˴���ʾ��ϸ˵��
%%ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
% R = abs(fun_SCMN(X));
%R_KA:����Э����
[M,N]=size(X);
rou_ba_t = 0;
for i = 1:N
    Ri = fun_Positive(X(:,i),5);
    rou_ba_t = rou_ba_t+norm(logm(Ri),'fro')^2/(N^2);
end
rou_ba = rou_ba_t-norm(R,'fro')^2/N;%R_KA
alpha0 = rou_ba/(rou_ba+norm(R-R_KA,'fro')^2);
alpha0 = max(min(1,alpha0),0);
% rou_ba = sum(diag(X'*X).^2)/N^2-sum(sum(abs(R).^2))/N;%��18��ʽ,
% alpha0 = rou_ba/(rou_ba+sum(sum(abs(R-R_KA).^2)));
R_CC = (1-alpha0)*R+alpha0*R_KA;
end


