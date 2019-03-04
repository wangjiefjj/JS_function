%%方阵的导向矢量
clc
clear 
close all
N_az = 12;
N_el = 12;
Az0 = 90;
El0 = 90;
d_az = 1.5;
d_el = 1.5;
lambda = 2*pi;
% sx = exp(1j*(0:N_az-1)*d_az*sin(El0/180*pi)*cos(Az0/180*pi));
% sy = exp(1j*(0:N_el-1)*d_el*sin(El0/180*pi)).';
% s = sy*sx;
% s = reshape(s,N_az*N_el,1);
s = fun_RectSteer(N_az, N_el, d_az, d_el, lambda,Az0, El0);
Az = 0:1:180;
El = 0:1:180;
A = zeros(length(El),length(Az));
for Eli = 1:length(El)
    for Azi = 1:length(Az)
        Gt = fun_Ax( N_az, N_el, d_az, d_el, lambda, Az(Azi), El(Eli), Az0, El0);
        Gr = fun_Ax( N_az, N_el, d_az, d_el, lambda, Az(Azi), El(Eli), Az0, El0);
        G = Gt*Gr/(N_az^2)/(N_el^2);        
%         sxi = exp(1j*(0:N_az-1)*d_az*sin(El(Eli)/180*pi)*cos(Az(Azi)/180*pi));
%         syi = exp(1j*(0:N_el-1)*d_el*sin(El(Eli)/180*pi)).';
%         si = G*reshape(syi*sxi,N_az*N_el,1)';
        si = G*fun_RectSteer(N_az, N_el, d_az, d_el, lambda,Az(Azi), El(Eli))';
        A(Eli,Azi)=si*s;
    end
end
[X,Y]=meshgrid(Az,El);
mesh(X,Y,(abs(A)))%db
xlabel('X方位角')
ylabel('Y俯仰角')