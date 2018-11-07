function [Mos] = fun_Mos(varargin)%%IC,Train,p,H,opt,(N,L,rho)
%FUN_IC 此处显示有关此函数的摘要
%   此处显示详细说明
%%
%%针对三种环境的模型选择
%%IC：字符选项，模型选择规则,AIC，AICc，ABIC，GIC
%%:Train：训练数据，
%%：p：目标导向适量
%%opt：适用数据选择，opt=1，只使用辅助数据，opt=2，适用辅助数据和主数据，
%%opt=3，只使用主数据 
%%H：字符选项，假设，H1,1假设均匀环境，H2,2假设部分均匀，H3,3假设，SIRP
%%N,数据维数
%%L，辅助数据个数

%%,H_num,:估计的参数个数
IC = varargin{1};
Train = varargin{2};
p = varargin{3};
H = varargin{4};
opt = varargin{5};
%% s函数计算 
switch(H)
    case 'H1'
        [s,H_num] = fun_s_H1(Train,p,opt);
    case 'H2' 
        [s,H_num] = fun_s_H2(Train,p,opt);
    case 'H3'
        [s,H_num] = fun_s_H3(Train,p,opt);     
end
s = -abs(s);
switch(IC)
    case 'AIC'
        Mos = fun_AIC(s,H_num); 
    case 'GIC'
       rho = varargin{6};
       Mos = fun_GIC(s,H_num,rho);
    case 'AICc'
       N = varargin{6};
       L = varargin{7};
       Mos = fun_AICc(s,H_num,N,L);
    case 'ABIC'
       L = varargin{6};
       Mos = fun_ABIC(s,H_num,L);   
end

end

