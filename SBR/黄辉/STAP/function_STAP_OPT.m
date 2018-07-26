function STAP_OPT= function_STAP_OPT(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,R_real,Ss,K,CNR)
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

 %% ��ѵ���������ƴ����ͨ���Ӳ�Э�������        
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
%% �õ�����ӦȨֵW_T      
    St = exp(j*[0:K-1].'*pi*fd(n));
    S = kron(St,Ss);   
    W_T = (inv(R)*S);
    miu=1/(S'*W_T);
%% �õ���������IF_OPT    
    SCNR_out = (abs(W_T'*S))^2/(W_T'*R_real*W_T);
    SCNR_in = (S'*S)/trace(R_real);%
    IF_OPT(n) = SCNR_out/SCNR_in; 
    W_OPT(:,n)=W_T; 
    Output(n)=abs(miu*W_T'*Data(:,Rangegate_now));
end
STAP_OPT.IF_OPT=IF_OPT;
STAP_OPT.W_OPT=W_OPT;
STAP_OPT.Output=Output;