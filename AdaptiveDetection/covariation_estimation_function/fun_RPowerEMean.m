function [ R ] = fun_RPowerEMean( X,alpha )
%FUN_RPOWEREMEAN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%Power-E�ļ��ξ�ֵЭ����
if nargin ==1
    alpha = 2;
end
[N,L] = size(X);
R_alpha = zeros(N,N);
for i = 1:L
     Ri = fun_Positive(X(:,i),4);
%      [UA, LA] = eig(Ri);
%      Ri = UA * (LA^alpha) * UA';
     Ri = Ri^alpha;
     R_alpha = R_alpha + Ri/L;
end
alpha = 1/alpha;
R = R_alpha^alpha;
% [UA, LA] = eig(R_alpha);
% R = UA * (LA^alpha) * UA';
end

