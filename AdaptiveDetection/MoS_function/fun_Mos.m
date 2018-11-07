function [Mos] = fun_Mos(varargin)%%IC,Train,p,H,opt,(N,L,rho)
%FUN_IC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%
%%������ֻ�����ģ��ѡ��
%%IC���ַ�ѡ�ģ��ѡ�����,AIC��AICc��ABIC��GIC
%%:Train��ѵ�����ݣ�
%%��p��Ŀ�굼������
%%opt����������ѡ��opt=1��ֻʹ�ø������ݣ�opt=2�����ø������ݺ������ݣ�
%%opt=3��ֻʹ�������� 
%%H���ַ�ѡ����裬H1,1������Ȼ�����H2,2���貿�־��ȣ�H3,3���裬SIRP
%%N,����ά��
%%L���������ݸ���

%%,H_num,:���ƵĲ�������
IC = varargin{1};
Train = varargin{2};
p = varargin{3};
H = varargin{4};
opt = varargin{5};
%% s�������� 
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

