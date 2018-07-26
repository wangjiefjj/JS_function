function STAP_OPT= function_STAP_OPT(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,R_real,Ss,K,CNR)
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

 %% 用训练样本估计待检测通道杂波协方差矩阵        
    TrainData= Data(:,TrainData_Range);
    R1 = TrainData*TrainData'/length(TrainData_Range);
    %     R = R1/max(max(abs(R1)))*CNR+eye(size(R1));
    %     R1=eye(192);
    %     R = CNR*R1./sum(eig(R1)/192);R = R+1*eye(192);
    %R = R1;    
    %R_real= R;
%%
 R=R_real ;

Nfd = length(fd);
for n = 1:Nfd
    %n
%% 得到自适应权值W_T      
    St = exp(j*[0:K-1].'*pi*fd(n));
    S = kron(St,Ss);   
    W_T = (inv(R)*S);
    miu=1/(S'*W_T);
%% 得到改善因子IF_OPT    
    SCNR_out = (abs(W_T'*S))^2/(W_T'*R_real*W_T);
    SCNR_in = (S'*S)/trace(R_real);%
    IF_OPT(n) = SCNR_out/SCNR_in; 
    W_OPT(:,n)=W_T; 
    Output(n)=abs(miu*W_T'*Data(:,Rangegate_now));
end
STAP_OPT.IF_OPT=IF_OPT;
STAP_OPT.W_OPT=W_OPT;
STAP_OPT.Output=Output;