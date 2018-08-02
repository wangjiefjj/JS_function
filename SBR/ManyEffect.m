%% 杂波带来的各种影响
clc
clear
close all

%%%距离模糊
%% 杂波协方差
isRu = 0;
isRotation = 0;
P00= fun_GenerateComplexTrainData(isRu,isRotation);
Rcn_00 = P00.R{1};
% [Rcn_00,El0,d,lambda] = fun_GenerateComplexR(isRu,isRotation);%无自转
% [Rcn_00,El0,d,lambda,fr,M,N] = fun_JWR(isRu,isRotation);%无自转
isRu = 0;
isRotation = 1;
P01= fun_GenerateComplexTrainData(isRu,isRotation);
Rcn_01 = P01.R{1};
% [Rcn_01] = fun_GenerateComplexR(isRu,isRotation);%有自转
% [Rcn_01] = fun_JWR(isRu,isRotation);%有自转
isRu = 1;
isRotation = 0;
P10= fun_GenerateComplexTrainData(isRu,isRotation);
Rcn_10 = P10.R{1};
% [Rcn_10] = fun_GenerateComplexR(isRu,isRotation);%有自转
% [Rcn_10] = fun_JWR(isRu,isRotation);%有自转
isRu = 1;
isRotation = 1;
P11= fun_GenerateComplexTrainData(isRu,isRotation);
Rcn_11 = P11.R{1};
Train11 = P11.TrainData;
N = length(Rcn_11);
Train11 = fun_TrainData('g',N,2*N, Rcn_11);
SCM11 = fun_SCMN(Train11);
% [Rcn_11] = fun_GenerateComplexR(isRu,isRotation);%有自转
% [Rcn_11] = fun_JWR(isRu,isRotation);%有自转
%%%%%%%%%%%%%%%%
%% 杂波子的秩
E=abs(eig(Rcn_00));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E,'r','LineWidth',2)
%%
E=abs(eig(Rcn_01));
E = 10*log10(sort(E,'descend')).';
plot(E,'k','LineWidth',2)
% % 
E=abs(eig(Rcn_10));
E = 10*log10(sort(E,'descend')).';
plot(E,'g','LineWidth',2)
% % 
E=abs(eig(Rcn_11));
E = 10*log10(sort(E,'descend')).';
plot(E,'b','LineWidth',2)
% %
E=abs(eig(SCM11));
E = 10*log10(sort(E,'descend')).';
plot(E,'y','LineWidth',2)
legend('无距离模糊无自转','无距离模糊有自转','有距离模糊无自转','有距离模糊有自转','有距离模糊有自转SCM')%
grid on
box on
xlabel('Index')
ylabel('Eig')
% axis([0,200,-1,60])
%% %% 杂波脊
% az = 0:1:180;     Laz = length(az);
% fd = -60928/2:60928/1000:60928/2;  Lfd = length(fd);
% fsp = d/lambda*sin(El0/180*pi)*cos(az*pi/180);
% omega = fd/fr;
% Pw2_00 = zeros(Lfd,Laz);
% for m=1:Laz
%     m
%     for n=1:Lfd
%         a = exp(1i*2*pi*fsp(m)*(0:N-1));                % Dummy Spatial Steering Vector.(Dummy虚拟)
%         b = exp(1i*2*pi*omega(n)*(0:M-1));              % Dummy Doppler Steering Vector
%         v = kron(b,a).';
%         Pw2_00(n,m) = abs(v'*Rcn_00*v)^2;                      %杂波脊
%         Pw2_01(n,m) = abs(v'*Rcn_01*v)^2;                      %杂波脊
%         Pw2_10(n,m) = abs(v'*Rcn_10*v)^2;                      %杂波脊
%         Pw2_11(n,m) = abs(v'*Rcn_11*v)^2;                      %杂波脊
%     end
% end
% %%Normalization:
% [Az,Doppler] = meshgrid(cos(az*pi/180),fd);
% max_value2 = max(max(Pw2_00));
% Pw2_00 = Pw2_00/max_value2;
% figure()
% colormap jet;
% imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw2_00)));
% shading interp;
% % xlim([-1 1])
% % ylim([-150 150]);
% str = ['cos(','\theta_{Az}',')'];
% xlabel(str);
% ylabel('Doppler Frequency (Hz)');
% h = colorbar;
% % h = colorbar('YTickLabel',{-80:10:0});
% set(get(h,'YLabel'),'String','Relative Power (dB)');
% title('无距离模糊无自转')
% % 
% % %%%
% % max_value2 = max(max(Pw2_01));
% Pw2_01 = Pw2_01/max_value2;
% figure()
% colormap jet;
% imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw2_01)));
% shading interp;
% % xlim([-1 1])
% % ylim([-150 150]);
% xlabel(str);
% ylabel('Doppler Frequency (Hz)');
% h = colorbar;
% % h = colorbar('YTickLabel',{-80:10:0});
% set(get(h,'YLabel'),'String','Relative Power (dB)');
% title('无距离模糊有自转')
% %%%
% % max_value2 = max(max(Pw2_10));
% Pw2_10 = Pw2_10/max_value2;
% figure()
% colormap jet;
% imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw2_10)));
% shading interp;
% % xlim([-1 1])
% % ylim([-150 150]);
% xlabel(str);
% ylabel('Doppler Frequency (Hz)');
% h = colorbar;
% % h = colorbar('YTickLabel',{-80:10:0});
% set(get(h,'YLabel'),'String','Relative Power (dB)');
% title('有距离模糊无自转')
% %%%
% % max_value2 = max(max(Pw2_11));
% Pw2_11 = Pw2_11/max_value2;
% figure()
% colormap jet;
% imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw2_11)));
% shading interp;
% % xlim([-1 1])
% % ylim([-150 150]);
% xlabel(str);
% ylabel('Doppler Frequency (Hz)');
% h = colorbar;
% % h = colorbar('YTickLabel',{-80:10:0});
% set(get(h,'YLabel'),'String','Relative Power (dB)');
% title('有距离模糊有自转')
%% MVD
fsp = 0;

% fsp = d/lambda/2 * sin(El0/180*pi);
a = exp(1i*2*pi*fsp*(0:N-1)).';    % Spatial Steering Vector.
L2 = 500;
omega = linspace(-0.5,0.5,L2);%% 归一化时间频率 
% fsp = d/lambda/2*cmj;
St=exp(1j*omega*(0:P.Np-1)).';
Ssr=exp(1j*fsp*(0:N-1)).';
Sst=exp(1j*fsp*(0:M-1)).';

iRcn_00 = inv(Rcn_00);
iRcn_01 = inv(Rcn_01);
iRcn_10 = inv(Rcn_10);
iRcn_11 = inv(Rcn_11);
for j = 1:L2
    St=exp(1j*omega*(0:P.Np-1)).';; % Temporal Steering Vector
    v = kron(Us.'*Rs.'*Us*Ut.'*Sst,kron(St,Ssr));;           % Space-TIme Steering Vector.
    MVD_00(j) = v'*iRcn_00*v;
    MVD_01(j) = v'*iRcn_01*v;
    MVD_10(j) = v'*iRcn_10*v;
    MVD_11(j) = v'*iRcn_11*v;
    
end
figure()
plot(omega*fr,10*log10(abs(MVD_00)/max(max(abs(MVD_00)))))
hold on 
plot(omega*fr,10*log10(abs(MVD_01)/max(max(abs(MVD_00)))),'r')
plot(omega*fr,10*log10(abs(MVD_10)/max(max(abs(MVD_00)))),'g')
plot(omega*fr,10*log10(abs(MVD_11)/max(max(abs(MVD_00)))),'k')
xlabel('多普勒(Hz)')
ylabel('相对功率(dB）')
legend('无距离模糊无自转','无距离模糊有自转','有距离模糊无自转','有距离模糊有自转')
grid on
%% %% FULLYSTAP
% isRu = 0;
% isRotation = 1;
% [Rcn_00,El0,d,lambda,fr,M,N] = fun_JWR(isRu,isRotation);
%% Target Space-Time Steering Vector
azt = 50; elt = 50;                                  % Target azimuth and elevation angles in degrees.
fdt = 20;                                               % Target Doppler Frequency.
omegact = fdt/fr;                                        % Normalized Target Frequency.
fspt = d/lambda*sin(elt*pi/180)*cos(azt*pi/180);     % Target Spatial Frequency.
at = exp(1i*2*pi*fspt*(0:N-1));                          % Target Spatial Steering Vector.
bt = exp(1i*2*pi*omegact*(0:M-1));                       % Target Doppler Steering Vector
vt = kron(bt,at).';                                      % Target Space-Time Steering Vector.

%% Optimum, Fully Adaptive STAP Solution
w00 = Rcn_00\vt;   
w01 = Rcn_01\vt;
w10 = Rcn_10\vt;
w11 = Rcn_11\vt;
%% Adapted Patterns
az = 0:.5:179;     Laz = length(az);
fd = -60928/2:60928/1000:60928/2;;  Lfd = length(fd);
fsp = d/lambda*sin(El0*pi/180)*cos(az*pi/180);
omega = fd/fr;
Pw1_00 = zeros(Lfd,Laz);
for m=1:Laz
    for n=1:Lfd
        a = exp(1i*2*pi*fsp(m)*(0:N-1));                % Dummy Spatial Steering Vector.(Dummy虚拟)
        b = exp(1i*2*pi*omega(n)*(0:M-1));              % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw1_00(n,m) = abs(w00'*v)^2;
        Pw1_01(n,m) = abs(w01'*v)^2;
        Pw1_10(n,m) = abs(w10'*v)^2;
        Pw1_11(n,m) = abs(w11'*v)^2;
        
    end
end

%% Normalization:
max_value00 = max(max(Pw1_00));
max_value01 = max(max(Pw1_01));
max_value10 = max(max(Pw1_10));
max_value11 = max(max(Pw1_11));
Pw00 = Pw1_00/max_value00;
Pw01 = Pw1_01/max_value01;
Pw10 = Pw1_10/max_value10;
Pw11 = Pw1_11/max_value11;


%% Plot the Adapted Pattern
figure();
[Az Doppler] = meshgrid(cos(az*pi/180),fd);
colormap jet;
% pcolor(Az, Doppler, 10*log10(abs(Pw)));
imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw00)));
hold on
plot(cos(azt*pi/180), fdt,'r'); 
view(0,-90)
shading interp;
% xlim([-1 1])
% ylim([-150 150]);
xlabel(str);
ylabel('多普勒 (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','相对功率 (dB)');
title('无距离模糊无自转')
%%%
figure();
colormap jet;
% pcolor(Az, Doppler, 10*log10(abs(Pw)));
imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw01)));
hold on
plot(cos(azt*pi/180), fdt,'r'); 
view(0,-90)
shading interp;
% xlim([-1 1])
% ylim([-150 150]);
xlabel(str);
ylabel('多普勒 (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','相对功率 (dB)');
title('无距离模糊有自转')
%%%
figure();
colormap jet;
% pcolor(Az, Doppler, 10*log10(abs(Pw)));
imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw10)));
hold on
plot(cos(azt*pi/180), fdt,'r'); 
view(0,-90)
shading interp;
% xlim([-1 1])
% ylim([-150 150]);
xlabel(str);
ylabel('多普勒 (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','相对功率 (dB)');
title('有距离模糊无自转')
%%%
figure();
colormap jet;
% pcolor(Az, Doppler, 10*log10(abs(Pw)));
imagesc(cos(az*pi/180), fd, 10*log10(abs(Pw11)));
hold on
plot(cos(azt*pi/180), fdt,'r'); 
view(0,-90)
shading interp;
% xlim([-1 1])
% ylim([-150 150]);
xlabel(str);
ylabel('多普勒 (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','相对功率 (dB)');
title('有距离模糊有自转')