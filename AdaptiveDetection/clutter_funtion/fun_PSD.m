function [ PSD,ft ] = fun_PSD( R,ft )
%PSD �˴���ʾ�йش˺�����ժҪ
% Capon PSD estimator
%%�� Estimation of the Covariance Matrix Based on Multiple A-Priori Models, A. De Maio��
%%N����һ�������յ���
%%R��Э�������
[M,~] = size(R);
if nargin == 1
    N = 1000;
    ft = linspace(-0.5,0.5,N+1);
elseif nargin == 2
    N = length(ft);
else
    error('��������һ��Э����')
end
nn = (0:M-1)';
PSD = zeros(N,1);
iR = inv(R);
for i = 1:N+1
    p = exp(1j*2*pi*nn*ft(i));
%     p = p/sqrt(M);
    PSD(i) = abs(p'*iR*p);
end
PSD = 10*log(PSD./max(abs(PSD)));
end

