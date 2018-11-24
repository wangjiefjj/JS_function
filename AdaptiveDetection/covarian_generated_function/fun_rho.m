function [ R_rho ] = fun_rho( rho,N,mi,fd,sigma2)
%FUN_ROU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%�������Э��������
%%%rho����������
%%%N������ʸ��ά��
%%mi: �η�
%%fd:�Ӳ�������
%%sigma2: Э����Ĺ���
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


