function [s] = fun_s_H2(Data,p,opt)
%%H2���貿�־��ȣ�s����
%%Data:���ݣ�
%%p:��������
%%opt��1ʱ�����ݽ��и������ݣ�2ʱ������Ϊ�������ݶ���
% %%m��Ҫ���Ƶò������� 
if nargin<2
    opt=1;
end
if opt == 1
    X=Data;
else
    X = Data(:,1:end-1);
    x0 = Data(:,end);
end
[N,K]=size(X);
%%����Э���������
% tau = abs(trace(X'*X)/N/K);
% Sigma = 1/K/tau*(X*X');
Sigma = eye(N);
for i = 1:10
    tau = abs(trace(X'/Sigma*X)/N/K);
    Sigma = 1/K/tau*(X*X');
end
%%����s
if opt == 1 %%ֻ�ø�������
    s = -N*K*log(pi)-N*K*log(tau)-K*log(det(Sigma))-trace(X'/Sigma*X)/tau;
elseif opt == 2  %%����������
    a = (p'/Sigma*x0)/(p'/Sigma*p);
%     tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
    tau0=1;
    s = -N*(K+1)*log(pi)-N*K*log(tau)-N*log(tau0)...   
        -(K+1)*log(det(Sigma))-trace(X'/Sigma*X)/tau...            
        -(x0-a*p)'/Sigma*(x0-a*p)/tau0;
elseif opt == 3 %%ֻ��������
    a = (p'/Sigma*x0)/(p'/Sigma*p);
%     tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
    tau0=1;
    s = -N*log(pi)-N*log(tau0)...   
        -log(det(Sigma))-(x0-a*p)'/Sigma*(x0-a*p)/tau0;   
end
end

