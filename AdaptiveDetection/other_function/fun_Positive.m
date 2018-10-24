function [ Rx ] = fun_Positive( X,opt )
%FUN_POSITIVE 此处显示有关此函数的摘要
%   此处显示详细说明
% 计算hermite正定矩阵（HPD)，7,8不是正定是结构话协方差
if nargin == 1
    opt = 1;
end
[N,~] = size(X);
% 
if opt==0
    Rx = X*X';
elseif opt == 1
%%《Covariance matrix estimation via geometric
%%barycenters and its application to radar training data selection》  
    X = X/sqrt(X'*X/N);%%归一化
    Rx = X*X';
    [V,D] = svd(Rx);
    Evalue = (diag(D));
    index_1 = find(Evalue<1);
    Evalue(index_1) = 1;
    % Evalue = sort(Evalue,'descend');
    D = diag(Evalue);
    Rx = abs(V)*D*abs(V');
elseif opt == 2
%%《Geometric barycenters for covariance estimation in compound-Gaussian
%%clutter》
    X = X/sqrt(X'*X/N);%%归一化
    Rx = X*X';
    [V,D] = svd(Rx);
    Evalue = (diag(D));
    Km=N;
%     Km = max(Evalue)/min(Evalue);
    xk = norm(X,'fro');
    lambdak = max(1, xk^2*Km/(Km^2+N-1) );
    Evalue(1)=Km*lambdak;
    Evalue(2:end)=lambdak;
%     Evalue(1) = lambdak * Km;
    D = diag(Evalue);
    Rx = abs(V)*D*abs(V');
%     Rx = Rx./xk^2;
elseif opt == 3
   X = X/sqrt(X'*X/N);%%归一化 
   Rx = (X*X'+ eye(N)); 
elseif opt == 4
   X = X/sqrt(X'*X/N);%%归一化
   Rx = fun_corrcoef(X) ;
elseif opt == 5
   Rx = fun_SFP(X,1);
elseif opt == 6
   Rx = fun_SFP(X,2);
elseif opt == 7  %%%persymmetric
    J = zeros(N,N);
    for i = 1:N
        J(i,N-i+1) = 1;
    end
    S = fun_NSCMN(X);
    Rx = 0.5*(S + J*conj(S)*J);  
elseif opt == 8 %%symmetric
    S = fun_NSCMN(X);
    Rx = real(S);  
end
end

