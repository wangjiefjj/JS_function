clc
clear 
close all
%%%%%%¹À¼Æ²ÎÊý
% Read_Display_Data
Data_process
load(matFile) 
% sig = sig.';
r = sig(:,Range);
[lambda,mu] = fun_IG_ML(r)
MTD = abs(fftshift(fft(sig,[],1)));
mesh(MTD)