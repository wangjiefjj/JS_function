function [ R ] = fun_RLogEMean( X, opt)
%%%Log-E的均值协方差
[N,L] = size(X);
logm_R = zeros(N,N);
% logm_Rr = zeros(N,N);
% logm_Ri = zeros(N,N);
if nargin<2
    opt = 4; %%正则化选项
end
for i = 1:L
    Ri = fun_Positive(X(:,i),opt);
%     if opt == 4 || opt == 1
       logm_R = logm_R + logm(Ri);
%     else
%        logm_R = logm_R + fun_Logm(Ri);
%     end  
% end
% if opt == 4 || opt == 1
    R = expm(logm_R/L);
% else
%     R = fun_Expm(logm_R/L);
% end


end

