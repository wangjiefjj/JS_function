clc
clear 
close all
%2017.09.15
% 
%%
Na=2;     % é˜µå…ƒæ•?
Np=4;     % è„‰å†²æ•?
N=Na*Np;
L=round(2*N); 
SNRout = -5:1:25; % è¾“å‡ºSNR
cos2=1;%%%Ê§Åä
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
% vt=randn(N,1)+1i*randn(N,1);
theta_sig=0.15;
nn=0:N-1;
vt=exp(-1i*2*pi*nn*theta_sig)'; %%%%%% µ¼ÏòÊ¸Á¿
rhoR=0.9;
fc = 0;%%ÔÓ²¨¶àÆÕÀÕ
for i=1:N
    for j=1:N
        R(i,j)=(rhoR^abs(i-j))*exp(1j*2*pi*fc*(i-j));
    end
end
iR=inv(R);
R_half=R^0.5;
% [Ddescend,Index]=sort(diagD,'descend'); figure; plot(10*log10(abs(Ddescend)),'b-x')
[UU,SS,VV]=svd(iR*vt);
vt_v=UU(:,2); %µ¼ÏòÊ¸Á¿ºÍÔÓ²¨µÄÕý½»·ÖÁ¿vt^H*iR*vt_v==0
weight=linspace(0,1,300);
for i=1:length(weight)
    vt_tmpt=weight(i)*vt+(1-weight(i))*vt_v;
    cos2_tmpt(i)=abs(vt_tmpt'*iR*vt).^2/abs(vt_tmpt'*iR*vt_tmpt*vt'*iR*vt);
end
[Min, Index]=min(abs(cos2-cos2_tmpt));
Weight=weight(Index);
vt_real=Weight*vt+(1-Weight)*vt_v;
figure; plot(weight,cos2_tmpt);
%% ÃÅÏÞ
Tamf = zeros(1,MonteCarloPfa);
Tglrt = zeros(1,MonteCarloPfa);
Tace = zeros(1,MonteCarloPfa);
Tabort = zeros(1,MonteCarloPfa);
Twabort = zeros(1,MonteCarloPfa);
Tprao = zeros(1,MonteCarloPfa);
Tdnamf = zeros(1,MonteCarloPfa);
Taed = zeros(1,MonteCarloPfa);
tic    
% h = waitbar(0,'Please wait...');
parfor i=1:MonteCarloPfa
    warning off
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
    X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
    R_half = R^(0.5);
    S=(R_half*X)*(R_half*X)';
    iS = inv(S);
    x = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    x = R_half*x;
    %%
    Tamf(i)=abs(vt'*iS*x)^2/abs(vt'*iS*vt);     %%%%%% AMFæˆ–è?wald
    tmp=abs(x'*iS*x);
    Tglrt(i)=Tamf(i)/(1+tmp);                   %%%%%% KGLRT
    Tace(i)=Tamf(i)/tmp;                        %%%%%% ACE
    Tabort(i)=(1+Tamf(i))/(2+tmp);              %%%%%% ABORT  % eq.(16) æ£?µ‹ç»Ÿè®¡é‡?
    Twabort(i)=1/(1+tmp)/(1-Tglrt(i))^2;        %%%%%% ABORT  % è§ä¼šè®®è®ºæ–‡ä¸­çš„eq.(18)
    Tace_bar=Tace(i)/(1-Tace(i));
    Tprao(i)=Tglrt(i)^2/(Tamf(i)*(1-Tglrt(i))); %%%%%% DMRao
    Tdnamf(i)=Tace_bar/tmp;                     %%%%%% DNAMF  % eq.(24) æ£?µ‹ç»Ÿè®¡é‡?
    Taed(i)=tmp;                                %%%%%% èƒ½é‡æ£?µ‹å™?

end
TACE=sort(Tace,'descend');
TKGLRT=sort(Tglrt,'descend');
TABORT=sort(Tabort,'descend');
TWABORT=sort(Twabort,'descend');
TDMRao=sort(Tprao,'descend');
TDNAMF=sort(Tdnamf,'descend');
TAED=sort(Taed,'descend');

Th_ACE=(TACE(floor(MonteCarloPfa*PFA-1))+TACE(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRT=(TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_ABORT=(TABORT(floor(MonteCarloPfa*PFA-1))+TABORT(floor(MonteCarloPfa*PFA)))/2;
Th_WABORT=(TWABORT(floor(MonteCarloPfa*PFA-1))+TWABORT(floor(MonteCarloPfa*PFA)))/2;
Th_DMRao=(TDMRao(floor(MonteCarloPfa*PFA-1))+TDMRao(floor(MonteCarloPfa*PFA)))/2;
Th_DNAMF=(TDNAMF(floor(MonteCarloPfa*PFA-1))+TDNAMF(floor(MonteCarloPfa*PFA)))/2;
Th_AED=(TAED(floor(MonteCarloPfa*PFA-1))+TAED(floor(MonteCarloPfa*PFA)))/2;

toc
% alpha=sqrt(SNRnum/abs(vt'*invR*vt)); % SNR=|alpha|^2*s'*R^(-1)*s 
a=0;b=2*pi;
%% ¿ªÊ¼¼ì²â
tic
counter_ace=0;
counter_glrt=0;
counter_abort=0;
counter_wabort=0;
counter_prao=0;
counter_dnamf=0;
counter_aed=0;
%%
alpha=sqrt(SNRnum/abs(vt'*iR*vt)); % SNR=|alpha|^2*s'*R^(-1)*s 
% alpha = sqrt(SNRnum/2);
%% ¼ì²â¸ÅÂÊ 
Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_ACE_mc = zeros(1,length(SNRout));
Pd_ABORT_mc = zeros(1,length(SNRout));
Pd_WABORT_mc = zeros(1,length(SNRout));
Pd_DMRao_mc = zeros(1,length(SNRout));
Pd_DNAMF_mc = zeros(1,length(SNRout));
Pd_AED_mc = zeros(1,length(SNRout));
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
        warning off
         X = (randn(N,L)+1i*randn(N,L))/sqrt(2);
         R_half = R^(0.5);
         S=(R_half*X)*(R_half*X)';
         iS=inv(S);
        x = (randn(N,1)+1i*randn(N,1))/sqrt(2);
        x = R_half*x;
        x = alpha(m)*vt+x;
       %%
        Tamf=abs(vt'*iS*x)^2/abs(vt'*iS*vt);  %%%%%% AMF
        tmp=abs(x'*iS*x);
        Tglrt=Tamf/(1+tmp);                   %%%%%% KGLRT
        Tace=Tamf/tmp;                        %%%%%% ACE
        Tabort=(1+Tamf)/(2+tmp);              %%%%%% ABORT  % eq.(16) æ£?µ‹ç»Ÿè®¡é‡?
        Twabort=1/(1+tmp)/(1-Tglrt)^2;        %%%%%% ABORT  % è§ä¼šè®®è®ºæ–‡ä¸­çš„eq.(18)
        Tace_bar=Tace/(1-Tace);
        Tprao=Tglrt^2/(Tamf*(1-Tglrt));       %%%%%% DMRao
        Tdnamf=Tace_bar/tmp;                  %%%%%% DNAMF  % eq.(24) æ£?µ‹ç»Ÿè®¡é‡?
        Taed=tmp;                             %%%%%% èƒ½é‡æ£?µ‹å™? 
        %%      
        if Tglrt>Th_KGLRT;      counter_glrt=counter_glrt+1;        end                
        if Tace>Th_ACE;         counter_ace=counter_ace+1;          end          
        if Tabort>Th_ABORT;     counter_abort=counter_abort+1;      end            
        if Twabort>Th_WABORT;   counter_wabort=counter_wabort+1;    end        
        if Tprao>Th_DMRao;      counter_prao=counter_prao+1;        end          
        if Tdnamf>Th_DNAMF;     counter_dnamf=counter_dnamf+1;      end          
        if Taed>Th_AED;         counter_aed=counter_aed+1;          end
    end
    Pd_KGLRT_mc(m)=counter_glrt/MonteCarloPd;        counter_glrt=0;
    Pd_ACE_mc(m)=counter_ace/MonteCarloPd;          counter_ace=0;
    Pd_ABORT_mc(m)=counter_abort/MonteCarloPd;      counter_abort=0;
    Pd_WABORT_mc(m)=counter_wabort/MonteCarloPd;    counter_wabort=0;
    Pd_DMRao_mc(m)=counter_prao/MonteCarloPd;        counter_prao=0;
    Pd_DNAMF_mc(m)=counter_dnamf/MonteCarloPd;      counter_dnamf=0;
    Pd_AED_mc(m)=counter_aed/MonteCarloPd;          counter_aed=0;  
   
end
close(h)
toc

figure(2);
hold on
plot(SNRout,Pd_KGLRT_mc,'b-+','linewidth',2)
plot(SNRout,Pd_ACE_mc,'r-x','linewidth',2)
plot(SNRout,Pd_ABORT_mc,'c-*','linewidth',2)
plot(SNRout,Pd_WABORT_mc,'m-P','linewidth',2)
plot(SNRout,Pd_DMRao_mc,'r-o','linewidth',2)
plot(SNRout,Pd_AED_mc,'r-d','linewidth',2)
plot(SNRout,Pd_DNAMF_mc,'g-s','linewidth',2); 
legend('KGLRT','ACE','ABORT','WABORT','DMRao','AED','DNAMF')
% legend({'KGLRT','AMF/DMwald','DMRao'},'FontSize',20)
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
% axis([5 25 0 1])
% clear TAMF TACE TKGLRT TABORT TWABORT TDMRao TDNAMF   Tamf Tace Tglrt Tabort Twabort Tprao Tdnamf Taed TAED X
