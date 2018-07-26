function [clutter_amp_fluct]=Fun_clutter_amp_fluct(clutter_amp_fluct_flag,reflect_cell_num,ampfluct_parameter,signal_parameter)
%%%%%%%%%%%%%%%%%%%杂波的幅度起伏模型%%%%%%%%%%%%%%%%%%
%1．Log-normal分布:
%                   y=1/(sigma*x*sqrt(2*pi))*exp(-(lnx-miu)^2/(2*sigma^2)),x>=0                   
%                   C_qf=lognrnd(miu,sigma,m,n);
%2．Weibull分布:
%                   y=a*b*x^(b-1)exp(-a*x^b), x>=0
%                   C_qf=wlbrnd(a,b,m,n); a,b分别为尺度参数和形状参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clutter_amp_fluct_flag==1
    
    model=ampfluct_parameter.model;
    %pulse_num=signal_parameter.pulse_num;
    %C_qf=ones(1,pulse_num);
    if model==1
        miu=ampfluct_parameter.miu;
        sigma=ampfluct_parameter.sigma;
        clutter_amp_fluct=lognrnd(miu,sigma,1,reflect_cell_num);
    end
    if model==2
        a=ampfluct_parameter.b;
        b=ampfluct_parameter.c;
        clutter_amp_fluct=wblrnd(a,b,1,reflect_cell_num);
    end
else
    clutter_amp_fluct=ones(1,reflect_cell_num);
end