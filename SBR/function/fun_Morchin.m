function [ sigma0 ] = fun_Morchin( graze, fo, opt, ss)
%FUN_MOR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%Morchinģ�ͣ�ģ����Ӳ�
%<<���Ӳ������л���Ԥ���״����þ������>>
% graze: �ӵؽ�
% fo: ��Ƶ
% opt: ����ѡ��
% ss: ����ȼ� 0~9
c = 3e8;
lambda = c/fo;
mu = sqrt(fo/10^9)/4.7;
if nargin == 3
    ss = 1;
end
if opt == 1 %%ɳĮ
    A = 0.00126;
    B = pi/2;
    beta0 = 0.14;
    phic = asin(lambda/4/pi/9.3/(beta0^2.2));
    if graze<phic
       sigmac0 = graze/phic;
    else
       sigmac0 = 1;
    end
    
elseif opt == 2 %%ũ��
    A = 0.004;
    B = pi/2;
    beta0 = 0.2;
    sigmac0 = 1;
elseif opt == 3 %%����
    A = 0.0126;
    B = pi/2;
    beta0 = 0.4;
    sigmac0 = 1;
elseif opt == 4 %%��ɽ
    A = 0.04;
    B = 1.24;
    beta0 = 0.5;
    sigmac0 = 1;
elseif opt == 5 %%����
    %<��(��)�Ӳ�������ģ���о�,����ި,����Ծ>
    A = 4e-7*10^(0.6*(ss+1));
    B = pi/2;
    beta0 = 2.44*((ss+1)^(1.08))/57.29;
    he = 0.025+0.046*ss^(1.72);
    phic = asin(lambda/4/pi/he);
    k = 1.9;
    if graze<phic
        sigmac0 = (graze/phic)^k;
    else
        sigmac0 = 1;
    end
    mu = 1;
end 
sigma0 = A*sigmac0*sin(graze)/lambda + ...
    mu*cot(beta0)^2*exp(-tan(B-graze)^2/tan(beta0)^2);
end

