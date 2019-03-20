clc
clear 
close all
%%%%%%¹À¼Æ²ÎÊý
% Read_Display_Data
%19980223_170435_ANTSTEP.CDF range = 8;
%19980204_224024_ANTSTEP.CDF range = 17;
[sig,range,matFile] = fun_Data_process(8,'19980223_170435_ANTSTEP.CDF');
load(matFile) 
t = [21:29];
offset = 50000;%50000;
parfor i = 1:10000%nsweep
    r = sig(i+offset,t);
    [lambda(i),mu(i)] = fun_IG_ML(r);
end
lambda_m=mean(lambda)
mu_m=mean(mu)
% MTD = abs(fftshift(fft(sig,[],1),1));
% mesh(MTD)