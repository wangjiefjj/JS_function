function [AIC] = fun_AIC(s,m)
%%AIC模型选择准则Akaike Information Criterion (AIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> 式（8）
%%s:s函数值
%%m：模型参数个数
AIC = -2*s+2*m;
end

