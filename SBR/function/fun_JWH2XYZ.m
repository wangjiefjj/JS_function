function [ X,Y,Z ] = fun_JWH2XYZ( longitude,latitude,h )
%%�ɵ�������ת��Ϊ���ĵع�����
longitude=longitude*pi/180;%����
latitude=latitude*pi/180;%γ��
% h %��
%����Ϊ��γ��+�����߶�
% a=6378.145;
% b=6356.76;
% e=0.08182;
% f=(a-b)/b;
% N=a/(1-e^2*(sin(longitude)^2));
% X=(N+h)*cos(longitude)*cos(latitude);
% Y=(N+h)*sin(longitude)*cos(latitude);
% Z=(N*a^2/b^2+h)*sin(latitude);
%����Ϊ�����׼��Բ�Ĳ���
%%%%������Բ����
X = h*cos(longitude)*cos(latitude);
Y = h*sin(longitude)*cos(latitude);
Z = h*sin(latitude);
end

