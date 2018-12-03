%% MOS��ʶ�𾫶�
clc
clear
close all
Data_process
% estiamtion
% matFile='19980204_224024_IPIX.mat';
load(matFile) 
Range = 19;
offset=0; %%ǰ10000������0����10000������50000
% lambda =  2.4072;
% mu = 1.3600;
%%%%��������
n = 2.5; %����������
%%%%�����������
SCNRout=0; % ���SNR
SCNRnum=10.^(SCNRout/10);
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
theta_sig = 0.2;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
a = sqrt(SCNRnum/abs(p'/R*p));
PFA=1e-3;% PFA=1e-4;
% MonteCarloPfa=1/PFA*100;
MC=nsweep;
L=round(n*N); 
Zhh = sig;
%%
H1_num = N^2+1;
H2_num = N^2+2;
H3_num = N^2+L+2;
tic
parfor i=1:10000-N+1
    %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
%     Train = fun_TrainData(str_train,N,L,Rc1,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
%     [x0,tau0] = fun_TrainData(str_train,N,1,Rc2,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
%     index_t1 = ceil(rand()*(M-10000))+2000;
    index_t1 = i+offset;
    Train1 = Zhh(index_t1:index_t1+N-1,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+N-1,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % �����źŽ������Ӳ�������    
    %%%����ź�
    x0=a*p+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
    %%%%%MOS%%%%%%%%%%%%
    s_H1 = -abs(fun_s_H1([Train,x0],p,2));
    s_H2 = -abs(fun_s_H2([Train,x0],p,2));
    s_H3 = -abs(fun_s_H3([Train,x0],p,2));
    %%%%ABIC%%%%%%%%%%%%%%
    H1_ABIC = fun_ABIC(s_H1,H1_num,L);
    H2_ABIC = fun_ABIC(s_H2,H2_num,L);
    H3_ABIC = fun_ABIC(s_H3,H3_num,L);
    Class_ABIC =[H1_ABIC,H2_ABIC,H3_ABIC];
    [~,Class_ABIC_num(i)] = min(Class_ABIC);
    %%%%GIC%%%%%%%%%%%%%%
    H1_GIC = fun_GIC(s_H1,H1_num,2);
    H2_GIC = fun_GIC(s_H2,H2_num,2);
    H3_GIC = fun_GIC(s_H3,H3_num,2);
    Class_GIC =[H1_GIC,H2_GIC,H3_GIC];
	[~,Class_GIC_num(i)] = min(Class_GIC);
	%%%%AIC%%%%%%%%%%%%%%
	H1_AIC = fun_AIC(s_H1,H1_num);
	H2_AIC = fun_AIC(s_H2,H2_num);
	H3_AIC = fun_AIC(s_H3,H3_num);
	Class_AIC =[H1_AIC,H2_AIC,H3_AIC];
	[~,Class_AIC_num(i)] = min(Class_AIC);  
end   
toc
Num = length(Class_AIC_num);
count_H1_ABIC = length(find(Class_ABIC_num==1))/Num*100;
count_H2_ABIC = length(find(Class_ABIC_num==2))/Num*100;
count_H3_ABIC = length(find(Class_ABIC_num==3))/Num*100;
count_H1_GIC = length(find(Class_GIC_num==1))/Num*100;
count_H2_GIC = length(find(Class_GIC_num==2))/Num*100;
count_H3_GIC = length(find(Class_GIC_num==3))/Num*100;
count_H1_AIC = length(find(Class_AIC_num==1))/Num*100;
count_H2_AIC = length(find(Class_AIC_num==2))/Num*100;
count_H3_AIC = length(find(Class_AIC_num==3))/Num*100;
y_H1 =[count_H1_AIC,count_H1_ABIC,count_H1_GIC];
y_H2 =[count_H2_AIC,count_H2_ABIC,count_H2_GIC];
y_H3 =[count_H3_AIC,count_H3_ABIC,count_H3_GIC];
y = [y_H1;y_H2;y_H3];
% y = y';
x={'H1','H2','H3'};
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','ABIC','GIC (\rho=2)');
% xlabel('Hy results')
ylabel('Percentage')
% set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gcf,'Position',[100 0 1200 1000])
grid on
box on
% str=['IPIX_recognition_Rang_',num2str(Range),'_',num2str(offset)];
% save(str);
