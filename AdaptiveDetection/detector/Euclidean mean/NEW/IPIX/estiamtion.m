clc
clear 
close all
%%%%%%���Ʋ���
% Read_Display_Data
Data_process
load(matFile) 
% sig = sig.';
r = sig(:,Range);
[lambda,mu] = fun_IG_ML(r)
MTD = abs(fftshift(fft(sig,[],1)));
mesh(MTD)