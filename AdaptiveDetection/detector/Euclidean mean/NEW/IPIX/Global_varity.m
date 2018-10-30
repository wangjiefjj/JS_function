clc
clear 
close all
%%%%%%%声明几个全局变量方便使用
cdfFile =  '19980205_185111_ANTSTEP.CDF';
Range = 14;
cdfFile_t = cdfFile;
if length(cdfFile_t)>17
    cdfFile_t(17:27)=[];
else
    cdfFile_t(4:7)=[];
end
matFile = [cdfFile_t,'IPIX.mat'];