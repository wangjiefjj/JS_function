clc
clear
close all
N = 8;
L = 40*N;
rho = 0.95;
R = fun_rho(rho,N,1);
% R = rand(N,N);
Train = fun_TrainData('g',N,L,R);
nn = 0:N-1;
theta_sig = 0.1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
s = s./sqrt(s'*s);% 归一化
K = N ; 
R_SCM = fun_SCMN(Train);
R_SCM_inv  = inv(R_SCM);
[R_SCM_LD, alpha] = fun_CC(Train,R_SCM,eye(N,N));
R_log = fun_RLogEMean(Train,2);
% R_log_LD = fun_LogCC_new(Train,eye(N,N));
r = 8;%rank(R_SCM);
T = zeros(N,r);
T(:,1) = s;
for i = 2:r
    t0 = R_SCM^(i-1)*s;
    while (1)
        p = t0;
        q = T'*p;
        v = T*q;
        t1 = p-v;
        if norm(t1)>norm(t0)/2
            break;
        else
            t0 = t1;
        end
    end
    t = t1/sqrt(t1'*t1);
    T(:,i) = t;
end
wcgs = T*inv(T'*R_SCM*T)*T'*s;
wscm = inv(R_SCM)*s;
wopt = inv(R)*s;
