function [s,m] = fun_s_H1(Data,p,opt)
%%H1假设均匀环境，s函数
%%Data:数据，
%%p:导向适量
%%opt：1时，数据仅有辅助数据，2时，数据为主副数据都有
% %%m：要估计得参数个数 
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
R = 1/K*(X*X');
% %% 如果有对称性的话
% J = zeros(N,N);
% for i = 1:N
%         J(i,N-i+1) = 1;
% end
% R = 0.5*(R + J*conj(R)*J); 
%%计算s
if opt == 1 %%只用辅助数据
    s = -N*K*log(pi)-K*log(det(R))-trace(X'/R*X);
    m = N^2;
elseif opt==2  %%主辅助数据
    a = (p'/R*x0)/(p'/R*p);
    s = -N*(K+1)*log(pi)-(K+1)*log(det(R))-trace(X'/R*X)-(x0-a*p)'/R*(x0-a*p);
    m = N^2+1;
elseif opt==3 %%只用主数据
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

