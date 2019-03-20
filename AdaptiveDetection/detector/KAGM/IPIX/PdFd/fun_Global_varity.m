%%%%%%%声明几个全局变量方便使用
function [matFile]=fun_Global_varity(str)
cdfFile =  str;
cdfFile_t = cdfFile;
cdfFile_t(17:27)=[];
matFile = [cdfFile_t,'IPIX.mat'];
end