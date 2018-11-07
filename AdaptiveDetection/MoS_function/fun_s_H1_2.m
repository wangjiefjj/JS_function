function [s] = fun_s_H1_2(Data,p,opt)
%%H1假设均匀环境，s函数
%%Data:数据，
%%p:导向适量
%%opt：1时，数据仅有辅助数据，2时，数据为主副数据都有
%%用Sigma 而不是R
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
%%计算协方差
Sigma = 1/K*(X*X');
for i = 1:10
    tau0 = abs(trace(X'/Sigma*X)/N/K);
    Sigma = 1/K/tau0*(X*X');
end
% R = R+eye(N);
%%计算s
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

