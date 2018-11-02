function [s] = fun_s_H2(Data,p,opt)
%%H2假设部分均匀，s函数
%%Data:数据，
%%p:导向适量
%%opt：1时，数据仅有辅助数据，2时，数据为主副数据都有
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
Sigma = 1/K*(X*X');
for i = 1:10
    tau = abs(trace(X'/Sigma*X)/N/K);
    Sigma = 1/K/tau*(X*X');
end
% Sigma = Sigma+(1/tau)*eye(N);
%%计算s
if opt == 1
    s = -N*K*log(pi)-N*K*log(tau)-K*log(det(Sigma))-trace(X'/Sigma*X)/tau;
elseif opt == 2
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0=1;
    s = -N*(K+1)*log(pi)-N*K*log(tau)-N*log(tau0)...   
        -(K+1)*log(det(Sigma))-trace(X'/Sigma*X)/tau...            
        -(x0-a*p)'/Sigma*(x0-a*p)/tau0;
elseif opt == 3
    a = (p'/Sigma*x0)/(p'/Sigma*p);
    tau0 = ((x0-a*p)'/Sigma*(x0-a*p))/N;
%     tau0=1;
    s = -N*log(pi)-N*log(tau0)...   
        -log(det(Sigma))-(x0-a*p)'/Sigma*(x0-a*p)/tau0;    
end
% real_s = real(s);
% if real_s>=0
%     s =abs(s);
% else
%     s =-abs(s);
% end
end

