function [ deltaR ] = fun_deltaR( H, R, Rs, Tr)
%FUN_DR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% H: ���Ǹ߶� km
% R: �ؾ� km
% Rs: б�� km
% dR: km
c = 3e8;  %����m/s
%% ���ģ������
Graze = fun_GrazeAngle(H,R,Rs)/180*pi;
deltaR =  c*Tr/2 *sec(Graze);
end

