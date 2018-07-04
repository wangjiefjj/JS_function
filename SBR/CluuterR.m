%% 按照JW的细节仿杂波, 线阵，无自转
clc
clear 
close all
% ha = 0; %%品台高度
% dphi = 0; %%波束宽度
%% 雷达系统参数
fo = 450e6; %9600e6 450e6     % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
M = 18;  %16                  % Number of Pulses per CPI.
Tp = 200e-6;                  % Pulse Width in sec.
N = 18;  %16                    % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 360;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d_hat = 1*lambda/2; %1   lambda/2         % Interelement Spacing
d = d_hat/(lambda/2);           % 归一化间隔
% Azimuth angle in degrees 方位角:
phi = linspace(-180,179,Nc);
Lphi = length(phi);
f = zeros(1,Lphi);
AF = ones(1,Lphi);           % Array Factor pre-allocation. 天线阵因子

% Platform Parameters:
ha = 9e3;%9e3 %700e3                    % Platform altitude in meters.
Vp = 49.9654;
if ha > 500e3
    Vp = fun_Vp(ha);          % 平台速度m/s
end
beta = Vp*2/fr/d_hat;         %%混叠系数;% beta parameter.
alpha1 = 45/180*pi;           %卫星纬度
eta = 90/180*pi;              %卫星倾角

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
Rcik = 130e3;                % (clutter) range of interest in meters. %%斜距
if ha >=500e3
    Rcik = 1300e3;                % (clutter) range of interest in meters. %%斜距
end
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Elevation (look-down angle) in rad. Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
phia = 0;                     % Velocity Misalignment angle in degrees.
dR = c/2/B;                   % Radar Range Resolution in meters.
%SBR
if ha >=500e3
    R = (fun_Rs2R(ha, Rcik));
    psi = fun_GrazeAngle(ha,R,Rcik)/180*pi;
    theta = pi/2 - fun_ELAngle(ha,R)/180*pi;
end


%% Clutter-to-Noise Ratio (CNR) Calculation 杂噪比估计(计算)
% Calculate the Voltage Element Pattern:
for i =1:Lphi
    if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
    else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
    end
end
figure('NumberTitle', 'off','Name', ...
       'Figure 9. The element voltage pattern. A -30-dB backlobe level is assumed.');
polardb(phi*pi/180,10*log10(abs(f)),-60,'g');

% Calculate the Array Factor (AF) (Voltage): 
% 这是线阵的
steering_angle = 0; % Angle of beam steering in degrees.%%波束的指向
for k=1:Lphi  
    AF(k) = sum(exp(-1i*pi*d*(0:N-1)*(sin(phi(k)*pi/180) ...
                  - sin(steering_angle*pi/180))));
end

% Calculate the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Clutter Patch RCS Calculation: %%反射因子乘以面积 
PatchArea = Rcik*dphi*dR*sec(psi); %杂波面积 ,dphi,应该是波束宽度吧
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

% Calculate the Clutter to Noise Ratio (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);
Ksic = sigma2*diag(ksi);

%% Clutter Covariance Matrix Computations 杂波协方差计算


% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
Ita = d_hat/lambda*cos(theta);
% Spatial frequency of the k-th clutter patch:
fspc = Ita*sin(phi*pi/180);
% Normalized Doppler Frequency of the k-th clutter patch:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);

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
    a(:,k) = exp(1i*2*pi*fspc(k)*(0:N-1));    % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));           % Space-TIme Steering Vector.
end

Rc = Vc*Ksic*Vc';                            % Eq. (64)

Rn = sigma2*eye(M*N);

Ru = Rc + Rn;
figure(3)
mesh(abs(Ru))
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

%% MVD
fsp = 0;
a = exp(1i*2*pi*fsp*(0:N-1)).';    % Spatial Steering Vector.
L2 = 200;
omega = linspace(-0.5,0.5,L2);%% 归一化时间频率 
iRu = inv(Ru);
for j = 1:L2
    b = exp(1i*2*pi*omega(j)*(0:M-1)).'; % Temporal Steering Vector
    v = kron(b,a);           % Space-TIme Steering Vector.
    MVD(j) = v'*iRu*v;
end
figure()
plot(10*log10(abs(MVD)/max(max(abs(MVD)))))