%% ����Mos���ʹ�ü����
%% ���Ȼ�����H1��ʹ��KGLRT
%% ���־��Ȼ���(H2) ʹ��KGLRT
%% SIRP����(H3)ʹ�� ANMF-AML
clc
clear
close all
%%%%��������
n = 1.5; %����������
%%%%�����������
SNRout=0:1:35; % ���SNR
SNRnum=10.^(SNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10);
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
theta_sig = 0.3;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
PFA=1e-3;% PFA=1e-4;
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;  
L=round(n*N); 
%% ���ּ��軷�������޼���
%% H1
str_train = 'g';
%%�Ӳ�Э����
fc=0;
sigmaf = 0.03; %%�Ӳ���չ��
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc1_H1 = CNRnum * toeplitz(rc);
Rc1_H1 = Rc1_H1+ eye(N) ;%
Rc2_H1 = Rc1_H1;
tic
parfor i = 1:MonteCarloPfa
%     warning('off')
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,Rc1_H1);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau0] = fun_TrainData(str_train,N,1,Rc2_H1); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% KGLRT
    Tkglrt(i) = fun_KGLRT(R_SCM,x0,p);
    %%%%% AMF
    Tamf(i) = fun_AMF(R_SCM,x0,p);
    %%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,p);
    %%%%% ANMF_AML
    Tanmf_AML(i) = fun_ANMF(R_AML,x0,p);
end
% close(h)
toc
TKGLRT=sort(Tkglrt,'descend');
TAMF=sort(Tamf,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_AML=sort(Tanmf_AML,'descend');

Th_KGLRT_H1 = (TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_AMF_H1 = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_SCM_H1 =(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML_H1 =(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;
%% H2
str_train = 'g';
%%�Ӳ�Э����
fc=0;
sigmaf = 0.03; %%�Ӳ���չ��
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc1_H2 = CNRnum * toeplitz(rc);
Rc1_H2 = Rc1_H2+ eye(N) ;%+ eye(N)
Rc2_H2 = 0.1*Rc1_H2;
parfor i = 1:MonteCarloPfa
    warning('off')
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,Rc1_H2);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau0] = fun_TrainData(str_train,N,1,Rc2_H2); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% KGLRT
    Tkglrt(i) = fun_KGLRT(R_SCM,x0,p);
    %%%%% AMF
    Tamf(i) = fun_AMF(R_SCM,x0,p);
    %%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,p);
    %%%%% ANMF_AML
    Tanmf_AML(i) = fun_ANMF(R_AML,x0,p);
end
% close(h)
toc
TKGLRT=sort(Tkglrt,'descend');
TAMF=sort(Tamf,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_AML=sort(Tanmf_AML,'descend');

Th_KGLRT_H2 = (TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_AMF_H2 = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM_H2 =(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML_H2 =(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;
%% H3
str_train = 'p';
 %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
%2ÿ����Ԫ����ֵһ����Ϊ���־��Ȼ���
opt_train = 1;  
lambda = 1;%%%ԽС�Ǹ�˹Խ����
mu = 1;
%%�Ӳ�Э����
sigmaf = 0.03; %%�Ӳ���չ��
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc1_H3 = CNRnum * toeplitz(rc);
Rc1_H3 = Rc1_H3+ eye(N) ;%+ eye(N)
Rc2_H3 = Rc1_H3;
parfor i = 1:MonteCarloPfa
%     warning('off')
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,Rc1_H3,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau0] = fun_TrainData(str_train,N,1,Rc2_H3,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% KGLRT
    Tkglrt(i) = fun_KGLRT(R_SCM,x0,p);
    %%%%% AMF
    Tamf(i) = fun_AMF(R_SCM,x0,p);
    %%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,p);
    %%%%% ANMF_AML
    Tanmf_AML(i) = fun_ANMF(R_AML,x0,p);
end
% close(h)
toc
TKGLRT=sort(Tkglrt,'descend');
TAMF=sort(Tamf,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_AML=sort(Tanmf_AML,'descend');

Th_KGLRT_H3 = (TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_AMF_H3 = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM_H3=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML_H3=(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%��⿪ʼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
counter_kglrt=0;
counter_amf=0;
counter_nscm=0;
counter_aml=0;
counter_mos=0;

H1_num = N^2+1;
H2_num = N^2+3;
H3_num = N^2+L+2;

Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_AMF_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_AML_mc = zeros(1,length(SNRout));
Pd_MOS_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    r = rand(MonteCarloPd);
    for i=1:MonteCarloPd
%         warning off
        if r(i) <0.3333
            Rc1=Rc1_H1;
            Rc2=Rc2_H1;
            str_train = 'g';
            flag_H(m,i)=1;
        elseif r(i) >=0.3333 && r(i) <0.6666
            Rc1=Rc1_H2;
            Rc2=Rc2_H2;
            str_train = 'g';
            flag_H(m,i)=2;
        elseif r(i) >=0.6666
            Rc1=Rc1_H3;
            Rc2=Rc2_H3;
            str_train = 'p';
            flag_H(m,i)=3;
        end
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,Rc1,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        [x0,tau0] = fun_TrainData(str_train,N,1,Rc2,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCM(Train));  
        R_NSCM = (fun_NSCMN(Train));  
        R_AML = fun_AML(Train);      
        %%%����ź�
        x0=alpha(m)*p+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% KGLRT
        Th_kglrt = fun_KGLRT(R_SCM,x0,p);
        %%%%% AMF
        Th_amf = fun_AMF(R_SCM,x0,p);
        %%%%% ANMF_SCM
        Th_scm = fun_ANMF(R_SCM,x0,p);
        %%%%% ANMF_AML
        Th_aml = fun_ANMF(R_AML,x0,p);
        %%%%%MOS%%%%%%%%%%%%
        s_H1 = -abs(fun_s_H1([Train,x0],p,2));
        s_H2 = -abs(fun_s_H2([Train,x0],p,2));
        s_H3 = -abs(fun_s_H3([Train,x0],p,2));
        H1_ABIC = fun_ABIC(s_H1,H1_num,L);
        H2_ABIC = fun_ABIC(s_H2,H2_num,L);
        H3_ABIC = fun_ABIC(s_H3,H3_num,L);
        Class_ABIC =[H1_ABIC,H2_ABIC,H3_ABIC];
        [~,Class_ABIC_num(m,i)] = min(Class_ABIC);
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if flag_H(m,i)==1
            if Th_kglrt>Th_KGLRT_H1;          counter_kglrt=counter_kglrt+1;        end      
            if Th_amf>Th_AMF_H1;          counter_amf=counter_amf+1;        end                
            if Th_scm>Th_SCM_H1;       counter_nscm=counter_nscm+1;    end
            if Th_aml>Th_AML_H1;       counter_aml=counter_aml+1;    end
        elseif flag_H(m,i)==2
            if Th_kglrt>Th_KGLRT_H2;          counter_kglrt=counter_kglrt+1;        end      
            if Th_amf>Th_AMF_H2;          counter_amf=counter_amf+1;        end                
            if Th_scm>Th_NSCM_H2;       counter_nscm=counter_nscm+1;    end
            if Th_aml>Th_AML_H2;       counter_aml=counter_aml+1;    end
        elseif flag_H(m,i)==3
            if Th_kglrt>Th_KGLRT_H3;          counter_kglrt=counter_kglrt+1;        end      
            if Th_amf>Th_AMF_H3;          counter_amf=counter_amf+1;        end                
            if Th_scm>Th_NSCM_H3;       counter_nscm=counter_nscm+1;    end
            if Th_aml>Th_AML_H3;       counter_aml=counter_aml+1;    end    
        end
        %%MOS�ж�
        if Class_ABIC_num(m,i) == 1
           if Th_kglrt>Th_KGLRT_H1;          counter_mos=counter_mos+1;        end 
        elseif Class_ABIC_num(m,i) == 2
           if Th_kglrt>Th_KGLRT_H2;          counter_mos=counter_mos+1;        end 
        elseif Class_ABIC_num(m,i) == 3
            if Th_aml>Th_AML_H3;       counter_mos=counter_mos+1;    end
        end
    end   
    Pd_KGLRT_mc(m)=counter_kglrt/MonteCarloPd;           counter_kglrt=0;
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;           counter_amf=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_AML_mc(m)=counter_aml/MonteCarloPd;        counter_aml=0;
    Pd_MOS_mc(m)=counter_mos/MonteCarloPd;        counter_mos=0;
end
tic
close(h)
figure(1);
hold on
plot(SNRout,Pd_KGLRT_mc,'r','linewidth',2)
plot(SNRout,Pd_AMF_mc,'g','linewidth',2)
plot(SNRout,Pd_NSCM_mc,'K-*','linewidth',2)
plot(SNRout,Pd_AML_mc,'k','linewidth',2)
plot(SNRout,Pd_MOS_mc,'b','linewidth',2)
h_leg = legend('KGLRT','AMF','ANMF-SCM','ANMF-AML','MOS');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
str=['MOS_PD_',num2str(L)];
save(str)