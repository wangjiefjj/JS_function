function [AICc] = fun_AICc(s,m,N,K)
%%AICcģ��ѡ��׼��corrected AIC
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> ʽ��16��
%%s:s����ֵ
%%m��ģ�Ͳ�������
t=(K+1)*N/((K+1)*N-m-1);
AICc = -2*s+2*m*t;
end
