function [SINR,ft] = fun_SINR(R,R_real)
%FUN_SINR �˴���ʾ�йش˺�����ժҪ
% Rank-Constrained Maximum Likelihood Estimation of 
% Structured Covariance Matrices(44)ʽ
[M,~] = size(R);
N=100;
ft = linspace(-0.5,0.5,N+1);
nn = (0:M-1)';
SINR = zeros(N,1);
iR = inv(R);
iR_real = inv(R_real);
for i = 1:N+1
    p = exp(1j*2*pi*nn*ft(i));
%     p = p/sqrt(M);
    SINR(i) = abs(p'*iR*p)^2/abs(p'*iR*R_real*iR*p)/abs(p'*iR_real*p);
end
SINR = 10*log(SINR);
SINR =SINR(:);
end

