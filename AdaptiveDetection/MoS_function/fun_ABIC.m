function [ABIC] = fun_ABIC(s,m,K)
%%ABIC模型选择准则Asymptotic Bayesian Information Criterion (ABIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> 式（18）
%%s:s函数值
%%m：模型参数个数
%%K:辅助数据长度
ABIC = -2*s+m*log(K);
end

