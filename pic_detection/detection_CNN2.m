% Example: Convolutional Neural Nets
% ---------------------
% ```matlab
clc
clear
close all
load randSNR.mat

% train_x = double(reshape(train_x',28,28,60000))/255;
train_x = trainData/max(max(max(trainData)));
% test_x = double(reshape(test_x',28,28,10000))/255;
test_x = trainData/max(max(max(testData)));
%%换成one-hot编码
for i = 1: length(trainLabels)
    if trainLabels(i) == 1
        train_y(:,i) = [1;0];%%有目标
    else
        train_y(:,i) = [0;1];%%无目标
    end
end
for i = 1: length(testLabels)
    if testLabels(i) == 1
        test_y(:,i) = [1;0];
    else
        test_y(:,i) = [0;1];
    end
end

%% ex1 Train a 6c-2s-12c-2s Convolutional neural network 
%will run 1 epoch in about 200 second and get around 11% error. 
%With 100 epochs you'll get around 1.2% error
rand('state',0)
cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};
cnn = cnnsetup(cnn, train_x, train_y);

opts.alpha = 0.5;
opts.batchsize = 100;
opts.numepochs = 100;
cnn = cnntrain(cnn, train_x, train_y, opts);
[er, bad] = cnntest(cnn, test_x, test_y);
%plot mean squared error
figure; plot(cnn.rL);
assert(er<0.12, 'Too big error');