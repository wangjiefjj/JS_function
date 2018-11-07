function [s] = fun_s_H1_2(Data,p,opt)
%%H1������Ȼ�����s����
%%Data:���ݣ�
%%p:��������
%%opt��1ʱ�����ݽ��и������ݣ�2ʱ������Ϊ�������ݶ���
%%��Sigma ������R
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
Sigma = 1/K*(X*X');
for i = 1:10
    tau0 = abs(trace(X'/Sigma*X)/N/K);
    Sigma = 1/K/tau0*(X*X');
end
% R = R+eye(N);
%%����s
if opt == 1
    s = -N*K*log(pi)-K*N*log(tau0)-K*log(det(Sigma))-trace(X'/Sigma*X);
else
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    s = -N*(K+1)*log(pi)-N*(K+1)*log(tau0)-(K+1)*log(det(Sigma))...
        -trace(X'/Sigma*X)-(x0-a*p)'/Sigma*(x0-a*p);
end
% real_s = real(s);
% if real_s>=0
%     s =abs(s);
% else
%     s =-abs(s);
% end
end

