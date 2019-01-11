function [F_t,Fr_fy]=Fun_trans_pattern(cita0,phi0,natenna_angle,crab_theta,cita,phi,lambda,natenna_parameter)
%%%%%%%%%%%%%%%%%%%%%阵列天线方向图%%%%%%%%%%%%%%%%%%%%%
%阵列天线采用：arry_row_num行，arry_col_num列,阵元间距arry_distance<=lambda/2,发射方向(cita0,phi0),列阵元误差delta_d_fy、行阵元误差delta_d_fw
%       方向维加权：I_fw=chebwin(arry_col_num,Jq_db)
%       俯仰维加权：I_fy=chebwin(arry_row_num,Jq_db)
%       发射方向图：F(cita,phi)=sum(sum(I_fw(m)*I_fy(n)*exp(j*2*pi*arry_distance/lambda((m-1+delta_d_fw（m))*(cos(cita)*cos(phi)-cos(cita0)*cos(phi0)+(n-1+delta_d_fy（n))*(sin(phi)-sin(phi0))))))
%       接收俯仰方向图：F_fy(phi)=sum(I_fy(n)*exp(j*2*pi*arry_distance/lambda((n-1+delta_d_fy（n))*(sin(phi)-sin(phi0)))))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arry_distance=natenna_parameter.arry_distance;
arry_col_num=natenna_parameter.arry_col_num;
Wq_azimuth_db=natenna_parameter.Wq_azimuth_db;
%delta_d_col=natenna_parameter.delta_d_col;
arry_row_num=natenna_parameter.arry_row_num;
Wq_pitch_db=natenna_parameter.Wq_pitch_db;
%delta_d_row=natenna_parameter.delta_d_row;
% arry_col_err=delta_d_col^2*randn(1,arry_col_num);
% arry_row_err=delta_d_row^2*randn(1,arry_row_num);
if Wq_azimuth_db==0
    I_azimuth=ones(arry_col_num,1);
else
    I_azimuth=chebwin(arry_col_num,Wq_azimuth_db);
end
if Wq_pitch_db==0
    I_pitch=ones(arry_row_num,1);
else
    I_pitch=chebwin(arry_row_num,Wq_pitch_db);
end
Lm=length(cita);
Ln=length(phi);
F_t=zeros(Ln,Lm);
Fr_fy=zeros(1,Ln);
I=I_azimuth*I_pitch';
m=0:arry_col_num-1;
n=0:arry_row_num-1;
% Dw=arry_distance*((m)+arry_cell_err);
Dw=arry_distance*(m);
Dy=arry_distance*(n);
% Fr_fy(k)=abs(F_fy*I_pitch);
for k=1:Ln
    F_fy=exp(j*2*pi/lambda*((sin(phi(k))-sin(phi0+crab_theta))*Dy));
    Fr_fy(k)=abs(F_fy*I_pitch);
    for g=1:Lm
        F_fw=exp(j*2*pi/lambda*((cos(cita(g))*cos(phi(k))-cos(cita0+natenna_angle)*cos(phi0+crab_theta))*Dw));
        F=I.*(F_fw'*F_fy);
        F_t(k,g)=abs(sum(sum(F)));
        %F_t(k,g)=sum(sum(F));
    end
    
end