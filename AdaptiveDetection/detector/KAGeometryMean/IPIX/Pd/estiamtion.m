clc
clear 
close all
%%%%%%¹À¼Æ²ÎÊý
% Read_Display_Data
%19980223_170435_ANTSTEP.CDF range = 8;
%19980204_224024_ANTSTEP.CDF range = 17;
[sig,range,matFile] = fun_Data_process(8,'19980223_170435_ANTSTEP.CDF');
load(matFile) 
r = sig(:,range);
[lambda,mu] = fun_IG_ML(r)
MTD = abs(fftshift(fft(sig,[],1)));
mesh(MTD)