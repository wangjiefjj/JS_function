%% 依据Mos结果使用检测器
%% 均匀环境（H1）使用KGLRT
%% 部分均匀环境(H2) 使用KGLRT
%% SIRP环境(H3)使用 ANMF-AML
clc
clear
close all
Data_process
estiamtion
% matFile='19980204_224024_IPIX.mat';
load(matFile) 
Range = 18;
% lambda =  2.4072;
% mu = 1.3600;
%%%%参数设置
n = 2; %几倍的样本
%%%%假设参数设置
SCNRout=0:1:35; % 输出SNR
SCNRnum=10.^(SCNRout/10);
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
theta_sig = 0.1;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
PFA=1e-3;% PFA=1e-4;
MonteCarloPd=1e4;  
% MonteCarloPfa=1/PFA*100;
MonteCarloPfa=M-MonteCarloPd;
L=round(n*N); 
Zhh = sig;
%% 门限计算
tic
parfor i = 1:MonteCarloPfa-N+1
%     warning('off')
%     %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%     Train = fun_TrainData(str_train,N,L,Rc1_H3,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%     [x0,tau0] = fun_TrainData(str_train,N,1,Rc2_H3,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%     index_t1 = ceil(rand()*(M-10000))+2000;
    index_t1 = i;
    Train1 = Zhh(index_t1:index_t1+N-1,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+N-1,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Th_KGLRT = (TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_AMF = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML=(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%检测开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
counter_kglrt=0;
counter_amf=0;
counter_nscm=0;
counter_aml=0;
counter_mos=0;

H1_num = N^2+1;
H2_num = N^2+3;
H3_num = N^2+L+2;

Pd_KGLRT_mc = zeros(1,length(SCNRout));
Pd_AMF_mc = zeros(1,length(SCNRout));
Pd_NSCM_mc = zeros(1,length(SCNRout));
Pd_AML_mc = zeros(1,length(SCNRout));
Pd_MOS_mc = zeros(1,length(SCNRout));
% alpha=sqrt(SCNRnum/abs(p'/R*p)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
alpha=sqrt(SCNRnum*mu/(lambda-1));    
h = waitbar(0,'Please wait...');
tic
for m=1:length(SCNRout)
    waitbar(m/length(SCNRout),h,sprintf([num2str(m/length(SCNRout)*100),'%%']));
    for i=1:MonteCarloPd
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%         Train = fun_TrainData(str_train,N,L,Rc1,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%         [x0,tau0] = fun_TrainData(str_train,N,1,Rc2,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%         index_t1 = ceil(rand()*(M-10000))+2000;
        index_t1 = i+ M-MonteCarloPd-N+1;
        Train1 = Zhh(index_t1:index_t1+N-1,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+N-1,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCM(Train));  
        R_NSCM = (fun_NSCMN(Train));  
        R_AML = fun_AML(Train);      
        %%%检测信号
        x0=alpha(m)*p+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Th_kglrt>Th_KGLRT;          counter_kglrt=counter_kglrt+1;        end      
        if Th_amf>Th_AMF;          counter_amf=counter_amf+1;        end                
        if Th_scm>Th_NSCM;       counter_nscm=counter_nscm+1;    end
        if Th_aml>Th_AML;       counter_aml=counter_aml+1;    end    
        %%MOS判断
        if Class_ABIC_num(m,i) == 1
           if Th_kglrt>Th_KGLRT;          counter_mos=counter_mos+1;        end 
        elseif Class_ABIC_num(m,i) == 2
           if Th_kglrt>Th_KGLRT;          counter_mos=counter_mos+1;        end 
        elseif Class_ABIC_num(m,i) == 3
            if Th_aml>Th_AML;       counter_mos=counter_mos+1;    end
        end
    end   
    Pd_KGLRT_mc(m)=counter_kglrt/MonteCarloPd;           counter_kglrt=0;
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;           counter_amf=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_AML_mc(m)=counter_aml/MonteCarloPd;        counter_aml=0;
    Pd_MOS_mc(m)=counter_mos/MonteCarloPd;        counter_mos=0;
end
toc
close(h)
figure();
hold on
plot(SCNRout,Pd_KGLRT_mc,'r','linewidth',2)
plot(SCNRout,Pd_AMF_mc,'g','linewidth',2)
plot(SCNRout,Pd_NSCM_mc,'K-*','linewidth',2)
plot(SCNRout,Pd_AML_mc,'k','linewidth',2)
plot(SCNRout,Pd_MOS_mc,'b','linewidth',2)
h_leg = legend('KGLRT','AMF','ANMF-SCM','ANMF-AML','MOS');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SCNRout),max(SCNRout),0,1])
grid on
str=['PD_',matFile(10:end-9),'_',num2str(L)];
save(str)