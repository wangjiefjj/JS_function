clc
clear 
close all
%%%%%%���Ʋ���
% Read_Display_Data
Data_process
load(matFile) 
range = 8;
% sig = sig.';
r = sig(:,range);
[alpha,beta] = fun_IG_ML(r)
MTD = abs(fftshift(fft(sig,[],1)));
mesh(MTD)