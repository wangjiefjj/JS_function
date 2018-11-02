function [s] = fun_s_H1(Data,p,opt)
%%H1假设均匀环境，s函数
%%Data:数据，
%%p:导向适量
%%opt：1时，数据仅有辅助数据，2时，数据为主副数据都有
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
% R = R+eye(N);
%%计算s
if opt == 1
    s = -N*K*log(pi)-K*log(det(R))-trace(X'/R*X);
elseif opt==2
    a = (p'/R*x0)/(p'/R*p);
    s = -N*(K+1)*log(pi)-(K+1)*log(det(R))-trace(X'/R*X)-(x0-a*p)'/R*(x0-a*p);
elseif opt==3 
    a = (p'/R*x0)/(p'/R*p);
    s = -N*log(pi)-log(det(R))-(x0-a*p)'/R*(x0-a*p);
end
% real_s = real(s);
% if real_s>=0
%     s =abs(s);
% else
%     s =-abs(s);
% end
end

