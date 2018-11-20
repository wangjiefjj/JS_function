function [s,m] = fun_s_H2(Data,p,opt)
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
J = zeros(N,N);
for i = 1:N
        J(i,N-i+1) = 1;
end
Sigma = 1/K*(X*X');
Sigma = 0.5*(Sigma + J*conj(Sigma)*J); 
for i = 1:10
    tau = abs(trace(X'/Sigma*X)/N/K);
    Sigma = 1/K/tau*(X*X');
    Sigma = 0.5*(Sigma + J*conj(Sigma)*J); 
end
% Sigma = Sigma+(1/tau)*eye(N);
%%����s
if opt == 1 %%ֻ�ø�������
    s = -N*K*log(pi)-N*K*log(tau)-K*log(det(Sigma))-trace(X'/Sigma*X)/tau;
    m = N^2+1;
elseif opt == 2  %%����������
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0=1;
    s = -N*(K+1)*log(pi)-N*K*log(tau)-N*log(tau0)...   
        -(K+1)*log(det(Sigma))-trace(X'/Sigma*X)/tau...            
        -(x0-a*p)'/Sigma*(x0-a*p)/tau0;
    m = N^2+3;
elseif opt == 3 %%ֻ��������
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0=1;
    s = -N*log(pi)-N*log(tau0)...   
        -log(det(Sigma))-(x0-a*p)'/Sigma*(x0-a*p)/tau0;   
    m = N^2+2;
end
m;
% real_s = real(s);
% if real_s>=0
%     s =abs(s);
% else
%     s =-abs(s);
% end
end

