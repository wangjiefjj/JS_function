function [ R_rho ] = fun_rho( rho,N,mi,fd,sigma2)
%FUN_ROU 此处显示有关此函数的摘要
%   此处显示详细说明
%%最基本的协方产生。
%%%rho，迟滞因子
%%%N，导向矢量维数
%%mi: 次方
%%fd:杂波多普勒
%%sigma2: 协方差的功率
R_rho = zeros(N,N);
if nargin<3 
    fd = 0;
    mi = 1;
    sigma2=1;
elseif nargin<4
    fd = 0;
    sigma2=1;
elseif nargin<5
    sigma2=1;    
end
L = length(rho);
jj = 1:N;
ii = (1:N).'; 
J = repmat(jj,N,1);
I = repmat(ii,1,N);
for l = 1:L
    R_rho(:,:,l)=sigma2*(rho(l).^abs(I-J).^mi).*exp(1j*2*pi*fd*(I-J));
end
% for i = 1:L
%     R_rho(:,:,i) = R_rho(:,:,i) + 10^(-2)*eye(N);
% end
end


