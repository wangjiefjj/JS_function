function [ R ] = fun_RLogEMean( X, opt)
%%%Log-E的均值协方差
[N,L] = size(X);
logm_R = zeros(N,N);
% logm_Rr = zeros(N,N);
% logm_Ri = zeros(N,N);
if nargin<2
    opt = 1;
end
if opt == 1
    for i = 1:L
        Ri = fun_Positive(X(:,i),4);
        logm_R = logm_R + logm(Ri)/L;
    end
else
    for i = 1:L%%归一化的
        Ri = fun_Positive(X(:,i),5);
        logm_R = logm_R + logm(Ri)/L;
%         Rir = real(Ri);
%         Rii = imag(Ri);
%         logm_Rr = logm_Rr + logm(Rir)/L;
%         logm_Ri = logm_Ri + logm(Rii)/L;
    end
end

% % logm_R = logm_R/L;
% % R = fun_Expm(logm_R);
R = expm(logm_R);
% R = expm(logm_Rr) + 1j*expm(logm_Ri);
end

