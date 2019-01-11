function [Sigma]=Fun_clutter_reflect_model(r,phi,airborne_parameter,radar_parameter,signal_parameter,entironment_parameter)
%%%%%%%%%%%%%%%%%%%�Ӳ���Ԫ����ģ��%%%%%%%%%%%%%%%%%%
%�Ӳ���Ԫ�����delta_s(r)=r*cita_3db*delta_r/sqrt(2)  rΪб��
%���淴���ʣ�
%           ��٤��ģ��Sigma_0(cita_g)=Gama*sin(cita_g)+sigma_os*exp(-(pi/2-cita_g)^2/delta_sigma^2)
%                         Gama          ��ʾ���������йص�ϵ��
%                         cigma_os      Ϊ���淴��ϵ��
%                         cita_g        Ϊ���ؽ�
%                         delta_sigma   Ϊ���淴�������
%           Morchinģ��Sigma_0(cita_g)=A*sigma_c*sin(cita_g)/lambda+u*(ctan(belta_0))^2*exp(-(tan(B-cita_g))^2/(tan(belta_0))^2)
%                                                        A           B       belta_0       sigma_c       u                 kk                       
%                 scene_parameter=1              ���Ӳ�  F1         pi/2       F2          F3            1                 1.9 
%                 scene_parameter=2              ɳ  Į  0.000126   pi/2       0.14        F3            sqrt(f_z)/4.7     1   
%                 scene_parameter=3              ũ  ��  0.004      pi/2       0.2         1             sqrt(f_z)/4.7
%                 scene_parameter=4              ��  ��  0.0126     pi/2       0.4         1             sqrt(f_z)/4.7
%                 scene_parameter=5              ��  ɽ  0.04       1.24       0.5         1             sqrt(f_z)/4.7
%               F1=4*10^(0.6*(ss+1)-7��;
%               F2=2.44*(ss+1)^1.08/57.29;
%               F3=(cita_g/cita_c)^kk    cita_g<cita_c;
%               F3=1                     cita_c<cita_g;
%               cita_c=asin(lambda/(4*pi*h_e));
%               ���Ӳ���h_e=0.025+0.046*ss^1.72;
%               ���Ӳ���h_e=9.3*belta_0^2.2;
%               ssΪ����ȼ���1��2��3��4��5��
%               ���ؽǣ�sin(cita_g)=H/r-r/(2*re)��sin(cita_g) = H/r-(r^2-H^2)/(2*re*r); 
%               �����ǣ�sin(phi)=H/r+r/(2*re)  phi=asin(H/r+r/(2*re));
%�ء����Ӳ����״��Ч�������rCS��:Sigma=Sigma_0*delta_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���������Morchinģ��
H=airborne_parameter.H;
cita_3db=radar_parameter.cita_3db;
delta_r=signal_parameter.delta_r;
lambda=signal_parameter.lambda;
f_z=signal_parameter.fz;
scene_parameter=entironment_parameter.scene_parameter;
ss=entironment_parameter.ss;
re=entironment_parameter.re;
% cita_g=asin(H/r-r/(2*re));
cita_g=asin( H/r-(r^2-H^2)/(2*re*r));
% delta_s=r*cita_3db*delta_r/sqrt(2);
delta_s=r*cita_3db*delta_r/cos(phi);
if scene_parameter==1
    kk=1.9;
    h_e=0.025+0.046*ss^1.72;
    cita_c=asin(lambda/(4*pi*h_e));
    A=4*10^(0.6*(ss+1)-7);
    B=pi/2;
    belta_0=2.44*(ss+1)^1.08/57.29;    
    if cita_g<=cita_c
        sigma_c=(cita_g/cita_c)^kk;
    else
        sigma_c=1;
    end
    u=1;
end
if scene_parameter==2
    kk=1; belta_0=0.14;
    h_e=9.3*belta_0^2.2;
    cita_c=asin(lambda/(4*pi*h_e));
    A=0.000126;
    B=pi/2;       
    if cita_g<=cita_c
        sigma_c=(cita_g/cita_c)^kk;
    else
        sigma_c=1;
    end
    u=sqrt(f_z)/4.7;
end
if scene_parameter==3    
    A=0.004;
    B=pi/2;
    belta_0=0.2;    
    sigma_c=1;
    u=sqrt(f_z)/4.7; 
end
if scene_parameter==4    
    A=0.0126;
    B=pi/2;
    belta_0=0.4;    
    sigma_c=1;
    u=sqrt(f_z)/4.7; 
end
if scene_parameter==5    
    A=0.04;
    B=1.24;
    belta_0=0.5;    
    sigma_c=1;
    u=sqrt(f_z)/4.7; 
end
Sigma_0=A*sigma_c*sin(cita_g)/lambda+u*(cot(belta_0))^2*exp(-(tan(B-cita_g))^2/(tan(belta_0))^2);
Sigma=Sigma_0*delta_s;
