clc
clear 
close all
%%%%%%%��������ȫ�ֱ�������ʹ��
cdfFile =  '19980204_224024_ANTSTEP.CDF';
Range = 17;
cdfFile_t = cdfFile;
if length(cdfFile_t)>17
    cdfFile_t(17:27)=[];
else
    cdfFile_t(4:7)=[];
end
matFile = [cdfFile_t,'IPIX.mat'];