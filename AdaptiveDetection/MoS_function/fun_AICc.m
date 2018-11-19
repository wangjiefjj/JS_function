function [AICc] = fun_AICc(s,m,N,K)
%%AICc模型选择准则corrected AIC
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> 式（16）
%%s:s函数值
%%m：模型参数个数
t=(K+1)*N/((K+1)*N-m-1);
AICc = -2*s+2*m*t;
% AICc = fun_AIC(s,m)+2*(m+1)*(m+2)/(K-m-2);
end

