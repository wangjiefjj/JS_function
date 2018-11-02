function [GIC] = fun_GIC(s,m,rho)
%%GIC模型选择准则Generalized Information Criterion (GIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> 式（9）
%%s:s函数值
%%m：模型参数个数
%%rho: 超参数，一般2,4
GIC = -2*s+(1+rho)*m;
end

