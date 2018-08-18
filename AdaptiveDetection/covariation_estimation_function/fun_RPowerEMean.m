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
     Ri = fun_Positive(X(:,i),5);
     [UA, LA] = eig(Ri);
     Ri = UA * (LA^alpha) * UA';
%     if alpha == -1
%         Ri = inv(Ri);
%     else
%         Ri = Ri^alpha;
%     end 
     R_alpha = R_alpha + Ri/L;
end
alpha = 1/alpha;
% if alpha == -1
%    R = inv(R_alpha);
% else
%    R = R_alpha^alpha;
% end 

[UA, LA] = eig(R_alpha);
R = UA * (LA^alpha) * UA';
end

