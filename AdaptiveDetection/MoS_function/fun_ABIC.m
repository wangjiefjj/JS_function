function [ABIC] = fun_ABIC(s,m,K)
%%ABICģ��ѡ��׼��Asymptotic Bayesian Information Criterion (ABIC)
%%<Model Order Selection Rules for Covariance
%%Structure Classification in Radar> ʽ��18��
%%s:s����ֵ
%%m��ģ�Ͳ�������
%%K:�������ݳ���
ABIC = -2*s+m*log(K);
end

