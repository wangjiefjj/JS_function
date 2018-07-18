%% 按照JW的细节仿杂波, 线阵
clc
clear 
close all
isRotation=1;
% ha = 0; %%品台高度
% dphi = 0; %%波束宽度
%% 雷达系统参数
fo = 9600e6; %9600e6 450e6     % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 500;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
M = 32;  %16                  % Number of Pulses per CPI.
Tp = 200e-6;                  % Pulse Width in sec.
N = 32;  %16                    % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 359;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % 间隔
% Azimuth angle in degrees 方位角:
% Az = linspace(-0,179,Nc);
Az = 0:1:179;
Nc = length(Az);
LAz = length(Az);
f = zeros(1,LAz);
AF = ones(1,LAz);           % Array Factor pre-allocation. 天线阵因子

% Platform Parameters:
H = 700e3;%9e3 %700e3                    % Platform altitude in meters.
Vp = 49.9654;
if H > 500e3
    Vp = fun_Vp(H);          % 平台速度m/s
end
beta = Vp*2/fr/d;         %%混叠系数;% beta parameter.
alpha1 = 30;           %卫星纬度deg
eta = 70;              %卫星倾角deg

%% Thermal Noise Power Computations 终端噪声计算
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power in Watts.

%% Clutter Patch Geometry computations 杂波路径几何结构计算
Rs = 130e3;                % (clutter) range of interest in meters. %%斜距
if H >=500e3
    Rs = 1300e3;                % (clutter) range of interest in meters. %%斜距
end
R = fun_Rs2R(H,Rs);
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
graze = fun_GrazeAngle(H,R,Rs);             % Grazing angle at the clutter patch in rad (flat earth model).
El = fun_ELAngle(H,R);                  % Elevation (look-down angle) in rad. Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
dR = c/2/B;                   % Radar Range Resolution in meters.
%SBR
% if H >=500e3
%     R = (fun_Rs2R(H, Rs));
%     graze = fun_GrazeAngle(H,R,Rs)/180*pi;
%     El = pi/2 - fun_ELAngle(H,R)/180*pi;
% end


%% Clutter-to-Noise Ratio (CNR) Calculation 杂噪比估计(计算)
% Calculate the Voltage Element Pattern:
% for i =1:LAz
%     if (Az(i))>=0 && (Az(i))<180
%         f(i) = sin(Az(i)*pi/180);
%     else
%         f(i) = 10^(be/10)*sin(Az(i)*pi/180);
%     end
% end
% figure('NumberTitle', 'off','Name', ...
%        'Figure 9. The element voltage pattern. A -30-dB backlobe level is assumed.');
% polardb(Az*pi/180,10*log10(abs(f)),-60,'g');

% Calculate the Array Factor (AF) (Voltage): 
% 这是线阵的
steering_angle = 90; % Angle of beam steering in degrees.%%波束的指向
for k=1:LAz  
    AF(k) = sum(exp(-1i*pi*d*(0:N-1)*(cos(Az(k)*pi/180) ...
                  - cos(steering_angle*pi/180))).*sin(Az(k)*pi/180));
end
f = AF;
% Calculate the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Clutter Patch RCS Calculation: %%反射因子乘以面积 
PatchArea = Rs*dphi*dR*sec(graze/180*pi); %杂波面积 ,dphi,应该是波束宽度吧
sigma0 = gamma*sin(graze/180*pi);
sigma = sigma0*PatchArea;

% Calculate the Clutter to Noise Ratio (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rs^4);
index=find(ksi==max(ksi));
% ksi(1:index-1)=0;
% ksi(index+1:end)=0;
Ksic = sigma2*diag(ksi);

%% Clutter Covariance Matrix Computations 杂波协方差计算


% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
Ita = d/lambda*sin(El/180*pi);
% Spatial frequency of the k-th clutter patch:
fspc = Ita*cos(Az/180*pi);
% Normalized Doppler Frequency of the k-th clutter patch:
if isRotation==1
    CrabA = fun_CrabAngle(alpha1,eta,H);
    CrabM = fun_CrabMagnitude(alpha1,eta,H);
else
    CrabA = 0;
    CrabM = 1;
end
omegac = beta*Ita*CrabM*cos(Az*pi/180 + CrabA*pi/180);
mod(omegac(index),-0.5)*fr
mod(omegac(index),0.5)*fr
% if  ha >500e3
%     Ita = d_hat/lambda*sin(theta);
%     fspc = Ita*cos(phi*pi/180);
%     omegac = beta*Ita*cos(phi*pi/180 + phia*pi/180);
% end


% Clutter Steering Vector Pre-allocation:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fspc(k)*(0:N-1));   % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));           % Space-TIme Steering Vector.
end

Rc = Vc*Ksic*Vc';                            % Eq. (64)

Rn = sigma2*eye(M*N);

Rcn = Rc + Rn;
figure(3)
mesh(abs(Rcn))
%% 杂波的空间频率和时间频率图
% L1 = 100;
% L2 = 200;
% fsp = linspace(-0.5,0.5,L1);%% 归一化空间频率
% omega = linspace(-0.5,0.5,L2);%% 归一化时间频率 
% for i = 1:L1
%     i
%     for j = 1:L2
%         a = exp(1i*2*pi*fsp(i)*(0:N-1)).';    % Spatial Steering Vector.
%         b = exp(1i*2*pi*omega(j)*(0:M-1)).'; % Temporal Steering Vector
%         v = kron(b,a);           % Space-TIme Steering Vector.
%         P(i,j) = v'*Ru*v;
%     end
% end
% figure()
% [Y,X] = meshgrid(omega,fsp);
% mesh(X,Y,10*log10(abs(P)))
% xlabel('归一化多普勒频率 ');
% ylabel('归一化空间频率');
% zlabel('相对功率(dB)');
% figure()
% colormap jet;
% imagesc(fsp,omega,10*log10(abs(P)))
% xlabel(' 归一化空间频率');
% ylabel('归一化多普勒频率');
% view(0,-90)
% h = colorbar;
% set(get(h,'YLabel'),'String','Relative Power (dB)');


%%
E=abs(eig(Rcn));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E,'r')
%% MVD
fsp = 0;
% fsp = d/lambda/2 * sin(El0/180*pi);
a = exp(1i*2*pi*fsp*(0:N-1)).';    % Spatial Steering Vector.
L2 = 500;
omega = linspace(-0.5,0.5,L2);%% 归一化时间频率 
% fsp = d/lambda/2*cmj;
iRcn = inv(Rcn);
for j = 1:L2
    b = exp(1i*2*pi*omega(j)*(0:M-1)).'; % Temporal Steering Vector
    v = kron(b,a);           % Space-TIme Steering Vector.
    MVD(j) = v'*iRcn*v;
end
figure()
plot(omega*fr,10*log10(abs(MVD)/max(max(abs(MVD)))))