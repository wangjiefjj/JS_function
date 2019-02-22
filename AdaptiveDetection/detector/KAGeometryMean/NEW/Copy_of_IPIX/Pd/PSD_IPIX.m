clc
clear 
close all
% Read_Display_Data
str_IPIX = '19980223_170435_ANTSTEP.CDF';
str_IPIX_t = str_IPIX(1:16);
[sig,Range,matFile]=fun_Data_process(8,str_IPIX);
load(matFile)
load(matFile) 
%%%%��������
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
rou = 0.90;  %%Э����������ɵĳ�������
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 5;     % ������
N = Na*Np;
PFA=1e-3;% PFA=1e-4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
before = 80+N-1; %%ȥǰ��֡��Ϊ����Э����
Zhh = sig;

R_SCM = zeros(N,N);
R_CC = zeros(N,N);
R_ML = zeros(N,N);
R_NSCM = zeros(N,N);
R_ECC = zeros(N,N);
R_LogM = zeros(N,N);
R_LogCC = zeros(N,N);
R_P = zeros(N,N);
R_PCC = zeros(N,N);
MC = 100;
for i = 1:MC
    i
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10000))+2000;
    %%R_KA
    R_KA = zeros(N,N);
    for ii = 1:before-N+1
        x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
        R_KA = R_KA+x_tt*x_tt'/(before-N+1);
    end
    Train1 = Zhh(index_t1:index_t1+N-1,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+N-1,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % �����źŽ������Ӳ�������
    %%Э�������
%     R_SCM = R_SCM + fun_SCMN(Train)/MC;
    R_CC = R_CC + fun_CC(Train,fun_NSCMN(Train),R_KA)/MC;
    R_ML = R_ML + fun_MLalpha(Train,fun_NSCMN(Train),R_KA,x0)/MC;
    R_NSCM = R_NSCM + fun_NSCMN(Train)/MC;
    R_ECC = R_ECC + fun_CCIter2(Train,R_NSCM,R_KA)/MC;
    R_LogM = R_LogM + fun_RLogEMean(Train,4)/MC;
    R_LogCC = R_LogCC + fun_LogCC_new(Train,R_KA,4)/MC;
    R_P = R_P + fun_RPowerEMean(Train,-1,4)/MC;
    R_PCC = R_PCC + fun_PowerCC(Train,R_KA,-1,4)/MC;
end
% [PSD_SCM,ft] = fun_PSD(R_SCM);
[PSD_CC,ft] = fun_PSD(R_CC);
[PSD_ML] = fun_PSD(R_ML);
[PSD_NSCM] = fun_PSD(R_NSCM);
[PSD_ECC] = fun_PSD(R_ECC);
[PSD_LogM] = fun_PSD(R_LogM);
[PSD_LogCC] = fun_PSD(R_LogCC);
[PSD_P] = fun_PSD(R_P);
[PSD_PCC] = fun_PSD(R_PCC);
figure
hold on
plot(ft,PSD_CC,'c-');
plot(ft,PSD_ML,'r-.');
plot(ft,PSD_NSCM,'g-.');
plot(ft,PSD_ECC,'g-');
plot(ft,PSD_LogM,'b-');
plot(ft,PSD_LogCC,'b-.');
plot(ft,PSD_P,'m-');
plot(ft,PSD_PCC,'m-.');
grid on
h_leg = legend('CC','ML','NSCM','ECC','LogM','LogCC','P','PCC');