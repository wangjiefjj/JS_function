function [R_MESS,M] = fun_MESS(R)
%FUN_MESS �˴���ʾ�йش˺�����ժҪ
%%%������С��������ֵ����Э����ݡ�
%%<Number of Source Signal Estimation By Mean Squared Eigenvalue Error (MSEE)>
[V,D]=svd(R);
d = diag(D);
[N,~] = size(R);
MinZm=1e6;
M=0;
for m = 1:N-1
%     sigma_w = mean(d(fun_g(d,m,N)+1:end));
    sigma_w = mean(d(m+1:end));
    mZm = m*var(d(1:m))+m*sigma_w;
    if mZm<MinZm
        MinZm = mZm;
        M=m;
    end
end

% sigma_w = mean(d(fun_g(d,M,N)+1:end));
sigma_w = mean(d(M+1:end));
R_MESS = V(:,1:M)*diag(d(1:M))*V(:,1:M)'+sigma_w*eye(N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%�õ��ĺ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%����������ƽ
function v=fun_v(d,m,N)
    v = 0;
    for i = m:N 
        vsigma = mean(d(i:end)); %%ʽ��46��
        v = v+(d(i)-vsigma)^2; %%ʽ��45��
    end
    v = mean(v);%%ʽ��45��         
end
function g = fun_g(d,m,N)
    maxg=0;
    g=m;
    for i = m:N-1
        maxg_t = fun_v(d,i,N)/fun_v(d,i+1,N); %%ʽ47
        if maxg_t>maxg
            maxg=maxg_t;
            g=i;
        end
    end
end