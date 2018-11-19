function [s,m] = fun_s_H1(Data,p,opt)
%%H1������Ȼ�����s����
%%Data:���ݣ�
%%p:��������
%%opt��1ʱ�����ݽ��и������ݣ�2ʱ������Ϊ�������ݶ���
% %%m��Ҫ���Ƶò������� 
if nargin<3
    opt=1;
end
if opt == 1
    X=Data;
else
    X = Data(:,1:end-1);
    x0 = Data(:,end);
end
[N,K]=size(X);
%%����Э����
R = 1/K*(X*X');
% %% ����жԳ��ԵĻ�
% J = zeros(N,N);
% for i = 1:N
%         J(i,N-i+1) = 1;
% end
% R = 0.5*(R + J*conj(R)*J); 
%%����s
if opt == 1 %%ֻ�ø�������
    s = -N*K*log(pi)-K*log(det(R))-trace(X'/R*X);
    m = N^2;
elseif opt==2  %%����������
    a = (p'/R*x0)/(p'/R*p);
    s = -N*(K+1)*log(pi)-(K+1)*log(det(R))-trace(X'/R*X)-(x0-a*p)'/R*(x0-a*p);
    m = N^2+1;
elseif opt==3 %%ֻ��������
    a = (p'/R*x0)/(p'/R*p);
    s = -N*log(pi)-log(det(R))-(x0-a*p)'/R*(x0-a*p);
    m = N^2+1;
end
m;
% real_s = real(s);
% if real_s>=0
%     s =abs(s);
% else
%     s =-abs(s);
% end
end

