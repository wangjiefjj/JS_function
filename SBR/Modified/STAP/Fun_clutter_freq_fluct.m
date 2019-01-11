function [clutter_freq_fluct]=Fun_clutter_freq_fluct(clutter_freq_fluct_flag,reflect_cell_num,entironment_parameter,airborne_parameter,signal_parameter,radar_parameter)
%%%%%%%%%%%%%%%%%%%%%杂波的频谱分布%%%%%%%%%%%%%%%%%%%%
%雷达杂波模拟中所采用的许多频谱都是高斯型或渐近高斯型的
%   S_p(f)=1/(sqrt(2*pi)*sigma_f)*exp(-f^2/(2*sigma_f^2));   
%           sigma_f=2*sigma_v/lambda;
%                   sigma_v=sqrt(sigma_v1^2+sigma_v2^2);
%                       sigma_v1为风速引起的地、海杂波的速度起伏  
%                               海杂波速度起伏:sigma_v1=0.101*Vm
%                               地杂波速度起伏:sigma_v1=0.0066*Vm
%                       sigma_v2为雷达平台的运动带来的速度起伏
%                               sigma_v2=V*cita_3db/(2*sqrt(2*ln2))
%                               Vm为风速，V为雷达平台速度
%以delta_f为频谱采样间隔 ，对功率谱密度进行M次采样(M一般取5或7）:
%   delta_f=0.6*f_3db=1.2*sqrt(2*ln(2))*sigma_f
%   S_ps(k)=sqrt(delta_f*S_p((k-(M+1)/2)*delta_f)); k=1,2,...M,
%相关高斯时间序列:X(m)=sum(S_ps(k)*exp(j*2*pi((k-(M+1)/2)*delta_f)*m/fr));fd为载机运动产生的多普勒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    Vm=entironment_parameter.Vm;
    V=airborne_parameter.V;
    lambda=signal_parameter.lambda;
    cita_3db=radar_parameter.cita_3db;
    scene_parameter=entironment_parameter.scene_parameter;
    pulse_num=signal_parameter.pulse_num;
    fr=signal_parameter.fr;
    
if clutter_freq_fluct_flag==1    
    if scene_parameter==1
        sigma_v1=0.101*Vm;
    else
        sigma_v1=0.0066*Vm;
    end
    sigma_v2=V*cita_3db/(2*sqrt(2*log(2)));
    sigma_v=sqrt(sigma_v1^2+sigma_v2^2);
    sigma_f=2*sigma_v/lambda;
    delta_f=1.2*sqrt(2*log(2))*sigma_f;
    M=7;
    k=1:M;
    X=zeros(M,reflect_cell_num);
    f=(k-(M+1)/2)*delta_f;
    S_p=1/(sqrt(2*pi)*sigma_f)*exp(-f.^2/(2*sigma_f^2));
    S_ps=sqrt(delta_f*S_p);
    for nn=1:reflect_cell_num
        for mm=1:M
            RSTN=rand(1);
            A=2*RSTN-1;AA=A*A;
            RSTN=rand(1);
            B=2*RSTN-1;BB=B*B;
            while(AA+BB)>1
                RSTN=rand(1);
                A=2*RSTN-1;AA=A*A;
                RSTN=rand(1);
                B=2*RSTN-1;BB=B*B;
            end
            C=(AA-BB)/(AA+BB);
            S=2*A*B/(AA+BB);
            X(mm,nn)=(C+j*S)*S_ps(mm);
        end
    end
    
    m=1:pulse_num;
    %clutter_freq_fluct=ones(1,pulse_num);
    clutter_freq_fluct=X.'*exp(j*2*pi*f'*m/fr);
else
    clutter_freq_fluct=ones(reflect_cell_num,pulse_num);
end