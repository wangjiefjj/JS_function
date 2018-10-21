function [Tparglrt] = fun_Partial_GER_GLRT(Train,x0,H)
%FUN_PARTIAL_GER_GLRT 此处显示有关此函数的摘要
%   此处显示详细说明
%%部分均匀+GER下的GLRT
[N,L] = size(Train);
S = fun_NSCM(Train);
H = S^(-0.5)*H;
x0 = S^(-0.5)*x0;
theta = (H'*H)\(H'*x0);
x1 = x0 - H*theta;
x1 = S^(-0.5)*x1;
P_H = H*(H'*H)^(-1)*H';
P_pH = eye(N)-P_H;
r = x0'*P_pH*x0;
fai1 = x1'*x1;
fai0 = x0'*x0;
% %%%%%%%%%%%%%%
% alpha1=N/(L+1-N)/(x0'*P_pH*x0);
% alpha0=N/(L+1-N)/((x0')*x0);
% %%%%%%%%%
% %%%%%%%%%%%
% t0 = alpha0^(-N/(L+1))*(1+alpha0*(x0')*x0);%
% t1 = alpha1^(-N/(L+1))*(1+alpha1*x0'*P_pH*x0);%
% Tparglrt = (t0/t1);
Tparglrt=((x0')*x0)/(x0'*P_pH*x0);
end

