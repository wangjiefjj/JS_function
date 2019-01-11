function STAP_JDL= function_STAP_JDL(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,Ts,Ss,Wq,K,CNR)
%% 
        %Data����ǰ������
        %fd:������ͨ��
        %Rangegate_now:��ǰ�����������
        %Rangegate_end:��������
        %N_TrainData��ѵ��������
        %N_timechannel:ѡȡʱ��ͨ����
        %N_spacechannel:ѡȡ����ͨ����
        %R_real:����ⵥԪ��ʵ��R��
        %Ss:����ⵥԪ������ʸ��
        %Ts:����ͨ������ʸ������
        %Wq:ʱ��ͨ����Ȩ
        %K:������
%%  �õ�ѵ�������ľ�����TrainData_Range
NN = floor(N_TrainData/2);
if NN >= Rangegate_now  %�����ѵ����������ʱȡ�����ͨ���Ҳ�ѵ����������Ϊ���ѵ������
    TrainData_Range = [1:Rangegate_now-1,Rangegate_now+1:N_TrainData+1];
end
if Rangegate_now > NN&&Rangegate_now < Rangegate_end-NN 
TrainData_Range = [Rangegate_now-NN-1:Rangegate_now-1,Rangegate_now+1:Rangegate_now+N_TrainData-NN+1];
end
if Rangegate_end-Rangegate_now <= NN %���Ҳ�ѵ����������ʱȡ�����ͨ�����ѵ����������Ϊ�Ҳ�ѵ������
    TrainData_Range = [Rangegate_end-N_TrainData:Rangegate_now-1,Rangegate_now+1:Rangegate_end];
end
%% 
Nfd = length(fd);
for n = 1:Nfd
    n
%% �õ���ά����T
%     if n <= floor(N_timechannel/2) %����������ͨ������ʱȡ�����ͨ���Ҳ�ͨ������Ϊ���ͨ��
%         fdn = n+floor(N_timechannel/2)+(1:(ceil(N_timechannel/2)-n));
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(fdn) fd(1:n+floor(N_timechannel/2))]);
%     end    
%     if n > floor(N_timechannel/2)&&n <= Nfd-floor(N_timechannel/2) 
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(n-floor(N_timechannel/2):n+floor(N_timechannel/2))]);
%     end
%     if n > Nfd-floor(N_timechannel/2) %���Ҳ������ͨ������ʱȡ�����ͨ�����ͨ������Ϊ�Ҳ�ͨ��
%         fdn = n-floor(N_timechannel/2)-(1:floor(N_timechannel/2)-Nfd+n);
%         Tt = diag(Wq)*exp(j*pi*(0:K-1).'*[fd(n-floor(N_timechannel/2):end) fd(fdn)]);
%     end   
    Tt = diag(Wq)*exp(j*pi*(0:K-1).'*(fd(n)+(-floor(N_timechannel/2):1:floor(N_timechannel/2))/(K)));
    T = kron(Tt,Ts);%��ά���󣨾����ά���ϵ���ʸ����
 %% �õ���ά���R��  
    %��ѵ���������ƴ����ͨ���Ӳ�Э�������
    Data_T = T'*Data;
    TrainData= Data_T(:,TrainData_Range);
    R1 = TrainData*TrainData'/length(TrainData_Range);
    %R = R1/max(max(abs(R1)))*CNR+eye(size(R1));
%     R_real = CNR*R_real./sum(eig(R_real)/192);R_real = R_real+1*eye(192);
    R_real_T = T'*R_real*T;
    R = R1;%R=R_real_T; 
%% �õ���ά�������ӦȨֵW_T      
    St = exp(j*[0:K-1]'*pi*fd(n));
    S = kron(St,Ss);    
    S_T = T'*S;
    W_T = (inv(R)*S_T);
    miu=1/(S_T'*W_T);
%% �õ���������IF_3DT�����յ�����ӦȨֵW_3DT
    
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
