function [ lambda,mu ] = fun_IG_ML(r)
%%《Maximum Likelihood Estimation for Compound-Gaussian Clutter with Inverse Gamma Texture》
%%(15) (14) 
%%复合高斯invert Gamma 纹理参数ML估计
% lambda:形状参数,lambda
% mu：尺度参数,mu
Ns = length(r);
beta_t = 0.1:0.01:5; %%%IPIX
r2 = abs(r).^2;
L_beta_min = 1e10;
for i = 1:length(beta_t)
    t1 = sum(r2./(beta_t(i)*r2+1));
    t2 = sum(log(beta_t(i)*r2+1));
    L_beta = (Ns*beta_t(i)*t1)/(Ns-beta_t(i)*t1)-t2;
    L_beta = abs(L_beta);
    if L_beta<L_beta_min
        L_beta_min = L_beta;
        mu = beta_t(i);
    end
end
t3 = sum(r2./(mu*r2+1));
lambda = Ns/(mu*t3)-1;
end

