function [AIC] = fun_AIC(s,m)
%%AICģ��ѡ��׼��Akaike Information Criterion (AIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> ʽ��8��
%%s:s����ֵ
%%m��ģ�Ͳ�������
AIC = -2*s+2*m;
end

