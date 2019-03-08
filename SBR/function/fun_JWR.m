function [ Rcn,El0,d,lambda,fr,M,N ] = fun_JWR( isRu, isRotation,Rs)
%%%JW的杂波协方差，不存在距离模糊
if nargin <3
    Rs = 1300e3;
end
%% 雷达系统参数
fo = 6000e6; %9600e6 450e6     % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 1000;%60928;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
M = 16;  %16                  % Number of Pulses per CPI.
Tp = 200e-6;                  % Pulse Width in sec.
N = 16;  %16                  % Number of Array Antenna Elements
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
Sys_DOF = M*N;

% Platform Parameters:
H = 500e3;%9e3 %700e3                    % Platform altitude in meters.
Vp = 49.9654;
if H > 500e3
    Vp = fun_Vp(H);          % 平台速度m/s
end
beta = Vp*2/fr/d;         %%混叠系数;% beta parameter.
alpha1 = 0;           %卫星纬度deg
eta = 45;              %卫星倾角deg

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
% Rs = 130e3;                % (clutter) range of interest in meters. %%斜距
% if H >=500e3
%     Rs = Rs;                % (clutter) range of interest in meters. %%斜距
% %     Rs = 1300e3;
% end
R = fun_Rs2R(H,Rs);
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
graze = fun_GrazeAngle(H,R,Rs);             % Grazing angle at the clutter patch in rad (flat earth model).
El = fun_ELAngle(H,R);                  % Elevation (look-down angle) in rad. Flat earth assumption.
El0 =El;
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
for i =1:LAz
    if (Az(i))>=0 && (Az(i))<180
        f(i) = sin(Az(i)*pi/180);
    else
        f(i) = 10^(be/10)*sin(Az(i)*pi/180);
    end
end
% f=ones(size(f));
% figure('NumberTitle', 'off','Name', ...
%        'Figure 9. The element voltage pattern. A -30-dB backlobe level is assumed.');
% polardb(Az*pi/180,10*log10(abs(f)),-60,'g');

% Calculate the Array Factor (AF) (Voltage): 
% 这是线阵的
steering_angle = 90; % Angle of beam steering in degrees.%%波束的指向
for k=1:LAz  
    AF(k) = sum(exp(-1i*pi*(d/lambda*2)*(0:N-1)*(cos(Az(k)*pi/180) ...
                  - cos(steering_angle*pi/180))));
end
% AF=ones(size(AF));
% f = AF;
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
%% 距离模糊
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
Rsu=c/(2*fr);                               %最大无模糊斜距m
Nk1 = floor((Rsmax-Rs)/Rsu);                %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                    %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                    %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
if isRu == 0
    Rsk = Rs;
else
    Rsk = [Rsk2,Rsk1];%Rs;%                     %模糊斜距m
end

Rk = fun_Rs2R(H,Rsk);                       %模糊地距m
elk = fun_ELAngle(H,Rk);                    %各模糊距离俯仰角deg
% grazek =  fun_GrazeAngle(H,Rk,Rsk)/180*pi;         %各模糊距离掠射角rad
%%%%%%%%%
Ita = d/lambda*sin(elk/180*pi);
% Normalized Doppler Frequency of the k-th clutter patch:
if isRotation == 1
    CrabA = fun_CrabAngle(alpha1,eta,H);
    CrabM = fun_CrabMagnitude(alpha1,eta,H);
else
    CrabA = 0;
    CrabM = 1;
end
% Spatial frequency of the k-th clutter patch:
for i = 1:length(Ita)
    fspc(i,:) = Ita(i)*cos(Az/180*pi);
    omegac(i,:) = beta*Ita(i)*CrabM*cos(Az*pi/180 + CrabA*pi/180);
end
% 
% mod(omegac(index),-0.5)*fr
% mod(omegac(index),0.5)*fr
% if  ha >500e3
%     Ita = d_hat/lambda*sin(theta);
%     fspc = Ita*cos(phi*pi/180);
%     omegac = beta*Ita*cos(phi*pi/180 + phia*pi/180);
% end


% Clutter Steering Vector Pre-allocation:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);
for i = 1:length(elk)
    i
    Vc_t = zeros(M*N,Nc);
    for k=1:Nc
        a(:,k) = exp(1i*2*pi*fspc(i,k)*(0:N-1));   % Spatial Steering Vector.
        b(:,k) = exp(1i*2*pi*omegac(i,k)*(0:M-1)); % Temporal Steering Vector
        Vc_t(:,k) = kron(b(:,k),a(:,k));           % Space-TIme Steering Vector.
    end
    Vc = Vc + Vc_t;
end

Rc = Vc*Ksic*Vc';                            % Eq. (64)

Rn = sigma2*eye(M*N);
% CNR_dB = 60;
% CNR = 10^(CNR_dB/10);
% Rn = CNR*Rc./sum(eig(Rn)/Sys_DOF);
Rcn = Rc + Rn;



end

