function [GIC] = fun_GIC(s,m,rho)
%%GICģ��ѡ��׼��Generalized Information Criterion (GIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> ʽ��9��
%%s:s����ֵ
%%m��ģ�Ͳ�������
%%rho: ��������һ��2,4
GIC = -2*s+(1+rho)*m;
end

