clc
clear 
close all
% 导向适量失配对检测性能的影响
%%
Na=2;     
Np=4;     
N=Na*Np;
optc = 'g';
opt_train = 1;%%1:SIRP,2:部分均匀
L=round(2*N); 
PFA=1e-3;% PFA=1e-4;
%%各种比
SNRout = 35; % 输出信噪比SNR
CNRout = 5; %杂噪比
JNRout = 5; %干噪比
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%虚警率和MC次数
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%信号
fs=[0.3,0.3175];
theta = ones(length(fs),1);
nn=(0:N-1)';
H = exp(-1i*2*pi*nn*fs); %%%%%% 导向矢量
vt = H*theta;

%% 各种协方差
%%杂波协方差
fc = 0;
sigmaf = 0.03; %%杂波谱展宽
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc = CNRnum * toeplitz(rc);
%%干扰协方差
fj = -0.2;
jam = exp(-1i*2*pi*nn*fj);
Rj = JNRnum * jam*jam';
%%检测单元协方差test
Rt = eye(N) + Rc + Rj;
%%参考单元协方差secondary 
Rs = eye(N) + Rc;
iRs = inv(Rs);
%% 失配导向矢量设置
cos2=0:0.1:1;%%%失配
iRt = inv(Rt);
[UU,SS,VV]=svd(iRt*vt);
s_v=UU(:,2); %%%%%% 与vt在白化空间正交，即：H^H*iR*s_v==0
P_H = H/(H'*iRt*H)*H';
for j = 1:length(cos2)
    weight=linspace(0,1,300);
    for i=1:length(weight)
        s_tmpt=weight(i)*vt+(1-weight(i))*s_v;
        cos2_tmpt(i)=abs(s_tmpt'*iRt*vt).^2/abs(s_tmpt'*iRt*s_tmpt*vt'*iRt*vt);
%         cos2_tmpt(i)=abs(s_tmpt'*iR*P_H*iR*s_tmpt)/abs(s_tmpt'*iR*s_tmpt);
    end
    [Min, Index]=min(abs(cos2(j)-cos2_tmpt));
    Weight=weight(Index);
    s_real(:,j)=Weight*vt+(1-Weight)*s_v;
end

%% 门限
Tpglrt = zeros(1,MonteCarloPfa);
Tprao = zeros(1,MonteCarloPfa);
Tpwald = zeros(1,MonteCarloPfa);
Tsglrt = zeros(1,MonteCarloPfa);
Tsrao = zeros(1,MonteCarloPfa);
Tswald = zeros(1,MonteCarloPfa);
tic    
% h = waitbar(0,'Please wait...'); 
parfor i=1:MonteCarloPfa
    warning off
    %% 数据产生
    Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
    x = fun_TrainData(optc,N,1,Rt ,3,1,opt_train);
    %% 检测器
    Tpglrt(i)=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% KGLRT
    Tprao(i)=fun_P_GRE_Rao(Train,x,H);                     %%%%%% RAO
    Tpwald(i)=fun_P_GRE_Wald(Train,x,H);                   %%%%%% Wald
    Tsglrt(i)=fun_SGLRT(Train,x,H);                   %%%%%% KGLRT
    Tsrao(i)=fun_SRAO(Train,x,H);                     %%%%%% RAO
    Tswald(i)=fun_SWALD(Train,x,H);                   %%%%%% Wald
 
end
toc
TPGLRT=sort(Tpglrt,'descend');
TPRAO=sort(Tprao,'descend');
TPWALD=sort(Tpwald,'descend');

TSGLRT=sort(Tsglrt,'descend');
TSRAO=sort(Tsrao,'descend');
TSWALD=sort(Tswald,'descend');

Th_PGLRT = (TPGLRT(floor(MonteCarloPfa*PFA-1))+TPGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_PRAO = (TPRAO(floor(MonteCarloPfa*PFA-1))+TPRAO(floor(MonteCarloPfa*PFA)))/2;
Th_PWALD=(TPWALD(floor(MonteCarloPfa*PFA-1))+TPWALD(floor(MonteCarloPfa*PFA)))/2;

Th_SGLRT = (TSGLRT(floor(MonteCarloPfa*PFA-1))+TSGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_SRAO = (TSRAO(floor(MonteCarloPfa*PFA-1))+TSRAO(floor(MonteCarloPfa*PFA)))/2;
Th_SWALD=(TSWALD(floor(MonteCarloPfa*PFA-1))+TSWALD(floor(MonteCarloPfa*PFA)))/2;
toc
a=0;b=2*pi;
h = waitbar(0,'Please wait...');
%% 检测概率 
Pd_PGLRT_mc = zeros(length(SNRnum),length(cos2));
Pd_PRAO_mc = zeros(length(SNRnum),length(cos2));
Pd_PWALD_mc = zeros(length(SNRnum),length(cos2));

Pd_SGLRT_mc = zeros(length(SNRnum),length(cos2));
Pd_SRAO_mc = zeros(length(SNRnum),length(cos2));
Pd_SWALD_mc = zeros(length(SNRnum),length(cos2));
%% 开始检测
for i_SNR = 1:length(SNRnum)
    counter_pglrt=0;
    counter_prao=0;
    counter_pwald=0;

    counter_sglrt=0;
    counter_srao=0;
    counter_swald=0;
    for m=1:length(cos2)
        waitbar(((i_SNR-1)*length(cos2)+m)/length(cos2)/length(SNRnum),...
            h,sprintf([num2str(((i_SNR-1)*length(cos2)+m)/length(cos2)/length(SNRnum)*100),'%%']));
       %% 信噪比
        alpha=sqrt(SNRnum(i_SNR)*N/abs(vt'*vt));
%         alpha=sqrt(SNRnum(i_SNR)/abs(vt'*iRs*vt));
        parfor i=1:MonteCarloPd
            warning off
            %% 数据产生
            Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
            x = fun_TrainData(optc,N,1,Rt,3,1,opt_train);
            x=alpha*s_real(:,m)+x;%+pp;    %%%%%%%  閲嶈  %%%%%%%%%%%%%
            %% 检测器
            Tpglrt=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% KGLRT
            Tprao=fun_P_GRE_Rao(Train,x,H);                     %%%%%% RAO
            Tpwald=fun_P_GRE_Wald(Train,x,H);                   %%%%%% Wald

            Tsglrt=fun_SGLRT(Train,x,H);                   %%%%%% SGLRT
            Tsrao = fun_SRAO(Train,x,H);               %%%SRAO
            Tswald = fun_SWALD(Train,x,H);               %%%SWald
            %%
            if Tpglrt>Th_PGLRT;         counter_pglrt=counter_pglrt+1;          end 
            if Tprao>Th_PRAO;         counter_prao=counter_prao+1;          end 
            if Tpwald>Th_PWALD;         counter_pwald=counter_pwald+1;          end 

            if Tsglrt>Th_SGLRT;         counter_sglrt=counter_sglrt+1;          end 
            if Tsrao>Th_SRAO;         counter_srao=counter_srao+1;          end 
            if Tswald>Th_SWALD;         counter_swald=counter_swald+1;          end 
        end
        Pd_PGLRT_mc(i_SNR,m)=counter_pglrt/MonteCarloPd;          counter_pglrt=0;
        Pd_PRAO_mc(i_SNR,m)=counter_prao/MonteCarloPd;          counter_prao=0;
        Pd_PWALD_mc(i_SNR,m)=counter_pwald/MonteCarloPd;          counter_pwald=0;
        Pd_SGLRT_mc(i_SNR,m)=counter_sglrt/MonteCarloPd;          counter_sglrt=0;
        Pd_SRAO_mc(i_SNR,m)=counter_srao/MonteCarloPd;          counter_srao=0;
        Pd_SWALD_mc(i_SNR,m)=counter_swald/MonteCarloPd;          counter_swald=0;
    end
end

close(h)
if length(SNRnum)==1
    figure(2);
    hold on
    plot(cos2,Pd_PGLRT_mc,'k-o','linewidth',2,'MarkerSize',2)
    plot(cos2,Pd_PRAO_mc,'r-x','linewidth',2,'MarkerSize',10)
    plot(cos2,Pd_PWALD_mc,'b-x','linewidth',2,'MarkerSize',10)

    plot(cos2,Pd_SGLRT_mc,'k-.o','linewidth',2,'MarkerSize',2)
    plot(cos2,Pd_SRAO_mc,'r-.x','linewidth',2,'MarkerSize',10)
    plot(cos2,Pd_SWALD_mc,'b-.x','linewidth',2,'MarkerSize',10)
    legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD')
%     legend('P-SGLRT','P-SRAO','P-SWALD')
%     legend('SGLRT','SRAO','SWALD')
    xlabel('cos^2{\phi}','FontSize',20)
    ylabel('Pd','FontSize',20)
    set(gca,'FontSize',20)
    grid on
elseif length(cos2)==1
    figure(3);
    hold on
    plot(SNRout,Pd_PGLRT_mc,'k-o','linewidth',2,'MarkerSize',2)
    plot(SNRout,Pd_PRAO_mc,'r-x','linewidth',2,'MarkerSize',10)
    plot(SNRout,Pd_PWALD_mc,'b-x','linewidth',2,'MarkerSize',10)

    plot(SNRout,Pd_SGLRT_mc,'k-.o','linewidth',2,'MarkerSize',2)
    plot(SNRout,Pd_SRAO_mc,'r-.x','linewidth',2,'MarkerSize',10)
    plot(SNRout,Pd_SWALD_mc,'b-.x','linewidth',2,'MarkerSize',10)
    legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD')
%     legend('P-SGLRT','P-SRAO','P-SWALD')
%     legend('SGLRT','SRAO','SWALD')
    xlabel('SNR/dB','FontSize',20)
    ylabel('Pd','FontSize',20)
    set(gca,'FontSize',20)
    grid on
else
    [X,Y] = meshgrid(SNRout,cos2);
    figure(4);

    hold on
    Llist = 0.3:0.2:0.99;
    contour(X,Y,Pd_PGLRT_mc','k-','linewidth',1,'ShowText','on',...
        'LevelList',Llist)
    contour(X,Y,Pd_PRAO_mc','r-.','linewidth',1,'ShowText','on',...
        'LevelList',Llist)
    contour(X,Y,Pd_PWALD_mc','b-.','linewidth',1,'ShowText','on',...
        'LevelList',Llist)

    % contour(X,Y,Pd_SGLRT_mc','k-.o','linewidth',1)
    % contour(X,Y,Pd_SRAO_mc','r-.x','linewidth',1)
    % contour(X,Y,Pd_SWALD_mc','b-.x','linewidth',1)
    legend('P-SGLRT','P-SRAO','P-SWALD')
    % legend('SGLRT','SRAO','SWALD')
    xlabel('SNR','FontSize',20)
    ylabel('cos^2{\phi}','FontSize',20)
    % ylabel('Pd','FontSize',20)
    set(gca,'FontSize',20)
    %     axis([15,30,0,1])
    grid on
end
    

