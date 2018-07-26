function STAP_JDL= function_STAP_JDL(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,Ts,Ss,Wq,K,CNR)
%% 
        %Data：当前景数据
        %fd:多普勒通道
        %Rangegate_now:当前待处理距离门
        %Rangegate_end:最大距离门
        %N_TrainData：训练样本数
        %N_timechannel:选取时域通道数
        %N_spacechannel:选取空域通道数
        %R_real:待检测单元真实的R阵
        %Ss:待检测单元空域导向矢量
        %Ts:空域通道导向矢量矩阵
        %Wq:时域通道加权
        %K:脉冲数
%%  得到训练样本的距离门TrainData_Range
NN = floor(N_TrainData/2);
if NN >= Rangegate_now  %在左侧训练样本不够时取待检测通道右侧训练样本来补为左侧训练样本
    TrainData_Range = [1:Rangegate_now-1,Rangegate_now+1:N_TrainData+1];
end
if Rangegate_now > NN&&Rangegate_now < Rangegate_end-NN 
TrainData_Range = [Rangegate_now-NN-1:Rangegate_now-1,Rangegate_now+1:Rangegate_now+N_TrainData-NN+1];
end
if Rangegate_end-Rangegate_now <= NN %在右侧训练样本不够时取待检测通道左侧训练样本来补为右侧训练样本
    TrainData_Range = [Rangegate_end-N_TrainData:Rangegate_now-1,Rangegate_now+1:Rangegate_end];
end
%% 
Nfd = length(fd);
for n = 1:Nfd
    n
%% 得到降维矩阵T
%     if n <= floor(N_timechannel/2) %在左侧多普勒通道不够时取待检测通道右侧通道来补为左侧通道
%         fdn = n+floor(N_timechannel/2)+(1:(ceil(N_timechannel/2)-n));
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(fdn) fd(1:n+floor(N_timechannel/2))]);
%     end    
%     if n > floor(N_timechannel/2)&&n <= Nfd-floor(N_timechannel/2) 
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(n-floor(N_timechannel/2):n+floor(N_timechannel/2))]);
%     end
%     if n > Nfd-floor(N_timechannel/2) %在右侧多普勒通道不够时取待检测通道左侧通道来补为右侧通道
%         fdn = n-floor(N_timechannel/2)-(1:floor(N_timechannel/2)-Nfd+n);
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(n-floor(N_timechannel/2):end) fd(fdn)]);
%     end   
    Tt = diag(Wq)*exp(j*pi*(0:K-1).'*(fd(n)+(-floor(N_timechannel/2):1:floor(N_timechannel/2))/(K)));
    T = kron(Tt,Ts);%降维矩阵（局域二维联合导向矢量）
 %% 得到降维后的R阵  
    %用训练样本估计待检测通道杂波协方差矩阵
    Data_T = T'*Data;
    TrainData= Data_T(:,TrainData_Range);
    R1 = TrainData*TrainData'/length(TrainData_Range);
    %R = R1/max(max(abs(R1)))*CNR+eye(size(R1));
%     R_real = CNR*R_real./sum(eig(R_real)/192);R_real = R_real+1*eye(192);
    R_real_T = T'*R_real*T;
    R = R1;%R=R_real_T; 
%% 得到降维后的自适应权值W_T      
    St = exp(j*[0:K-1]'*pi*fd(n));
    S = kron(St,Ss);    
    S_T = T'*S;
    W_T = (inv(R)*S_T);
    miu=1/(S_T'*W_T);
%% 得到改善因子IF_3DT和最终的自适应权值W_3DT
    
%     R_real_T = R;
    SCNR_out = (abs(W_T'*S_T))^2/(W_T'*R_real_T*W_T);
     SCNR_in = (S'*S)/trace(R_real);
%    SCNR_in =1/CNR;
    IF_JDL(n) = SCNR_out/SCNR_in; 
    W_JDL(:,n)=miu*T*W_T; 
    Output(n)=abs(W_JDL(:,n)'*Data(:,Rangegate_now));
end
STAP_JDL.IF_JDL=IF_JDL;
STAP_JDL.W_JDL=W_JDL;
STAP_JDL.Output=Output;
