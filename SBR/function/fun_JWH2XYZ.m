function [ X,Y,Z ] = fun_JWH2XYZ( longitude,latitude,h )
%%由地理坐标转换为地心地固坐标
longitude=longitude*pi/180;%经度
latitude=latitude*pi/180;%纬度
% h %高
%以上为经纬度+距地面高度
% a=6378.145;
% b=6356.76;
% e=0.08182;
% f=(a-b)/b;
% N=a/(1-e^2*(sin(longitude)^2));
% X=(N+h)*cos(longitude)*cos(latitude);
% Y=(N+h)*sin(longitude)*cos(latitude);
% Z=(N*a^2/b^2+h)*sin(latitude);
%以上为地球标准椭圆的参数
%%%%以下是圆地球
X = h*cos(longitude)*cos(latitude);
Y = h*sin(longitude)*cos(latitude);
Z = h*sin(latitude);
end

