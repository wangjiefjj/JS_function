function [s,m] = fun_s_H2(Data,p,opt)
%%H2假设部分均匀，s函数
%%Data:数据，
%%p:导向适量
%%opt：1时，数据仅有辅助数据，2时，数据为主副数据都有
% %%m：要估计得参数个数 
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
%%计算协方差和纹理
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
%%计算s
if opt == 1 %%只用辅助数据
    s = -N*K*log(pi)-N*K*log(tau)-K*log(det(Sigma))-trace(X'/Sigma*X)/tau;
    m = N^2+1;
elseif opt == 2  %%主辅助数据
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0=1;
    s = -N*(K+1)*log(pi)-N*K*log(tau)-N*log(tau0)...   
        -(K+1)*log(det(Sigma))-trace(X'/Sigma*X)/tau...            
        -(x0-a*p)'/Sigma*(x0-a*p)/tau0;
    m = N^2+3;
elseif opt == 3 %%只用主数据
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

