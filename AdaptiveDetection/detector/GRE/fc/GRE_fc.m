clc
clear 
close all
% 
%%
Na=2;     
Np=4;     
N=Na*Np;
optc = 'g';
opt_train = 1;%%1:SIRP,2:���־���
L=round(2*N); 
cos2=1;%%%ʧ��
PFA=1e-3;% PFA=1e-4;
%%���ֱ�
SNRout = 30; % ��������SNR
CNRout = 15; %�����
JNRout = 15; %�����
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%�龯�ʺ�MC����
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%�ź�
fs=[0.3];%,0.3175
theta = ones(length(fs),1);
nn=(0:N-1)';
H = exp(-1i*2*pi*nn*fs)/sqrt(N); %%%%%% ����ʸ��
vt = H*theta;

% [UU,SS,VV]=svd(iR*H);
% vt_v=UU(:,end); %����ʸ�����Ӳ�����������vt^H*iR*vt_v==0,
% q = 2*vt_v;%%GRE�е�q
% R0 = R + q*q';
a=0;b=2*pi;
%% �����
alpha=sqrt(SNRnum*N/abs(vt'*vt));
% alpha=sqrt(SNRnum/abs(vt'/Rt*vt)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
%% ������ 
Pd_SGLRT_mc = zeros(1,length(JNRout));
Pd_SRAO_mc = zeros(1,length(JNRout));
Pd_SWALD_mc = zeros(1,length(JNRout));
h = waitbar(0,'Please wait...');
fc = -0.5:0.025:0.5;
tic
for m=1:length(fc)
    waitbar(m/length(fc),h,sprintf([num2str(m/length(fc)*100),'%%']));
    %% ����Э����
    %%�Ӳ�Э����
    sigmaf = 0.03; %%�Ӳ���չ��
    rc =  exp(-1i*2*pi*nn*fc(m)-2*(nn*pi*sigmaf).^2);
    Rc = CNRnum * toeplitz(rc);
    %%����Э����
    fj = -0.2;
    jam = exp(-1i*2*pi*nn*fj);
    Rj = JNRnum * jam*jam';
    %%��ⵥԪЭ����test
    Rt = eye(N) + Rc + Rj;
    %%�ο���ԪЭ����secondary 
    Rs = eye(N) + Rc;
    %% ����
    Tsglrt = zeros(1,MonteCarloPfa);
    Tsrao = zeros(1,MonteCarloPfa);
    Tswald = zeros(1,MonteCarloPfa);
    tic    
    % h = waitbar(0,'Please wait...'); 
    for i=1:MonteCarloPfa
        warning off
        %% ���ݲ���
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt ,3,1,opt_train);
        %% �����
        Tsglrt(i)=fun_SGLRT(Train,x,H);                   %%%%%% KGLRT
        Tsrao(i)=fun_SRAO(Train,x,H);                     %%%%%% RAO
        Tswald(i)=fun_SWALD(Train,x,H);                   %%%%%% Wald
    end
    toc
    TSGLRT=sort(Tsglrt,'descend');
    TSRAO=sort(Tsrao,'descend');
    TSWALD=sort(Tswald,'descend');

    Th_SGLRT = (TSGLRT(floor(MonteCarloPfa*PFA-1))+TSGLRT(floor(MonteCarloPfa*PFA)))/2;
    Th_SRAO = (TSRAO(floor(MonteCarloPfa*PFA-1))+TSRAO(floor(MonteCarloPfa*PFA)))/2;
    Th_SWALD=(TSWALD(floor(MonteCarloPfa*PFA-1))+TSWALD(floor(MonteCarloPfa*PFA)))/2;
    %% ��ʼ���
    counter_sglrt=0;
    counter_srao=0;
    counter_swald=0;
    parfor i=1:MonteCarloPd
        warning off
        %% ���ݲ���
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt,3,1,opt_train);
        x=alpha*vt+x;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
       %%
        
        Tsglrt=fun_SGLRT(Train,x,H);                   %%%%%% SGLRT
        Tsrao = fun_SRAO(Train,x,H);               %%%PRAO
        Tswald = fun_SWALD(Train,x,H);               %%%PWald
        %%
        if Tsglrt>Th_SGLRT;         counter_sglrt=counter_sglrt+1;          end 
        if Tsrao>Th_SRAO;         counter_srao=counter_srao+1;          end 
        if Tswald>Th_SWALD;         counter_swald=counter_swald+1;          end 
    end
    Pd_SGLRT_mc(m)=counter_sglrt/MonteCarloPd;          counter_sglrt=0;
    Pd_SRAO_mc(m)=counter_srao/MonteCarloPd;          counter_srao=0;
    Pd_SWALD_mc(m)=counter_swald/MonteCarloPd;          counter_swald=0;
end
close(h)
toc

figure(2);
hold on
plot(fc,Pd_SGLRT_mc,'k-o','linewidth',2,'MarkerSize',2)
plot(fc,Pd_SRAO_mc,'r-x','linewidth',2,'MarkerSize',10)
plot(fc,Pd_SWALD_mc,'b-x','linewidth',2,'MarkerSize',10)
legend('GLRT','RAO','WALD')
xlabel('fc','FontSize',20)
ylabel('PD','FontSize',20)
set(gca,'FontSize',20)
grid on
box on
save('fc');