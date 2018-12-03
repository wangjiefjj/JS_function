function [s] = fun_s_H3(Data,p,opt)
%%H3����SIRP��s����
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
Sigma = eye(N);
for i = 1:5
    S = 0;
    for k = 1:K  
        tauk(k) = abs(X(:,k)'/Sigma*X(:,k)/N);
        S = S+ X(:,k)*X(:,k)'/tauk(k)/K;   
    end
    Sigma = S;
end
%%����s
S = 0;
for k = 1:K
    S = S+X(:,k)'/Sigma*X(:,k)/tauk(k);
end
if opt == 1
    s = -N*K*log(pi)-sum(N*log(tauk))-K*log(det(Sigma))-S;
    m = N^2+K;
elseif opt==2
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = (x0-a*p)'/Sigma*(x0-a*p)/N;
%     tau0 = 1;
    s = -N*(K+1)*log(pi)-sum(N*log(tauk))-N*log(tau0)...   
        -(K+1)*log(det(Sigma))-S...            
        -(x0-a*p)'/Sigma*(x0-a*p)/tau0;
elseif opt == 3
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0 = 1;
    s = -N*log(pi)-N*log(tau0)...   
        -log(det(Sigma))-(x0-a*p)'/Sigma*(x0-a*p)/tau0;
end
end

