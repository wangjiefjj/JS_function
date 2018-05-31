function [ wd ] = fun_Wd_beta0( H,R,Az,beta)
%FUN_WD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �޵�����ת�Ķ�����Ƶ��
Re = 6373e3; %flat earth radius��m
t1 = Re*sin(R./Re).*cos(Az);
t2 = fun_R2Rs(H,R);
wd = beta *t1./t2;
end
