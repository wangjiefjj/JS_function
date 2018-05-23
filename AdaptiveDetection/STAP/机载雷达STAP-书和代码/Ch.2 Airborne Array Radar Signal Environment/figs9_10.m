%% Figures 9 and 10. CNR per Column (element) as a function of azimuth angle.

%%
% Coded by Ilias Konsoulas, 16 Dec. 2014. 
% Code provided for educational purposes only. All rights reserved.

clc;  clear;  close all;

%% Radar System Operational Parameters.
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
dc = 0.06;                    % Duty Factor or Duty Cycle 6%
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4*1e6;                   % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
% Number of Pulses per CPI: M = 18
Tp = 200e-6;                  % Pulse Width in sec.
N = 18;                       % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 360;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:179;
Lphi = length(phi);
f = zeros(1,Lphi);
AF = zeros(1,Lphi);           % Array Factor pre-allocation.

%% Platform Parameters.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.

%% Thermal Noise Power Computations.
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.        
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts 

%% Clutter Patch Geometry Computations. 
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters. 
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor. 
% bwdth = 2*asin(0.446*lambda/N/d);  % Antenna Azimuth 3-dB beamwidth in rad.

%% Calculate the Voltage Element Pattern.
for i =1:Lphi
     if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
     else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
     end
end

%% Plot in polar coordinates the magnitude of the element voltage gain.
figure('NumberTitle', 'off','Name', ...
       'Figure 9. The element voltage pattern. A -30-dB backlobe level is assumed.');
polardb(phi*pi/180,10*log10(abs(f)),-60,'g');

%% Calculate and Plot the Array Factor (AF) (Voltage).
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
     AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ... 
                   - sin(steering_angle*pi/180))).*cos(phi(k)*pi/180));
end

figure('NumberTitle', 'off','Name','The voltage Array Factor for N=18 elements.','Position',[1 1 1000 400]);
subplot(1,2,1);
polardb(phi*pi/180,10*log10(abs(AF)),-60,'r')
subplot(1,2,2);
plot(phi, 10*log10(abs(AF)));
grid on;
ylim([-30     15]);
xlim([-180 180]);

%% Calculate and Plot the Full Array Transmit Power Gain.
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate and Plot the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Total Reception Gain:
Grec = 10^(Gr/10)*grgain;

% figure('NumberTitle', 'off','Name','The Array Transmit and Element Receive Power Gain ');
% subplot(1,2,1);
% polardb(phi*pi/180,10*log10(abs(grgain)),-90)
% subplot(1,2,2);
% plot(phi, 10*log10(abs(Gtgain)));
% grid on;
% ylim([-60 50]);
% xlim([-180 180]);

%% Clutter Patch RCS Calculation.
PatchArea = Rc*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

%% Calculate and Plot the Clutter to Noise Ration (CNR) for each clutter patch.
ksi = Pt*Gtgain.*Grec*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rc^4);   % Eq. (58)

figure('NumberTitle', 'off','Name','Figure 10. Received CNR per column as a function of azimuth. ','Position',[1 1 650 500]);
plot(phi, 10*log10(abs(ksi)),'LineWidth',1.5);
grid on;
ylim([-80 40]);
xlim([-180 180]);
ylabel('CNR (dB)');
xlabel('Azimuth Angle (deg)');
title('CNR as a function of Azimuth angle');
