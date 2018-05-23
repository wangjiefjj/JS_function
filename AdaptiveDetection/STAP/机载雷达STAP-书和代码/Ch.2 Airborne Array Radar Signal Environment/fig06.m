%% Figure 6. The â=1 Clutter Ridge.

%%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters.
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B = 4e6;                      % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
% M = 18;                     % Number of Pulses per CPI: 
Tp = 200e-6;                  % Pulse Width in sec.
N = 18;                       % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 361;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing.

% Azimuth angle in degrees:
phi = -180:.5:180;
Lphi = length(phi);

%% Platform Parameters:
va = 50;                      % Platform velocity in m/sec.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Depression angle to ik-th clutter patch (flat earth model).

%% Spatial Frequency of ik-th Clutter Patch:
fsp = d/lambda*cos(theta)*sin(phi*pi/180);         % Eq. (54)

% Doppler Frequency of ik-th Clutter Patch:
fd = 2*va/lambda*cos(theta)*sin(phi*pi/180);       % Eq. (69)

% Normalized Doppler Frequency:                    
omegac = 2*va*Tr/d*fsp;                            % Eq. (70)

beta = 2*va*Tr/d;                                  % Eq. (71)

%% Plot Doppler Frequency vs Azimuth and Normalized Doppler Frequency vs Spatial Frequency.
figure('NumberTitle', 'off','Name',...
    ['Figure 6. The â =1 clutter ridge for Side Looking Airborne Radar (SLAR). The PRF is ', ...
                                                        num2str(fr),' Hz.'],'Position', [50  50  900 400] );
subplot(1,2,1);
plot(sin(phi*pi/180),fd,'.');
grid on;
ylabel('Doppler Frequency f_c(\theta_c,\phi_c) (Hz)');
xlabel('sin(\phi_c)');

subplot(1,2,2);
plot(fsp,omegac,'.');
grid on;
ylabel('Normalized Doppler Frequency \omega_c');
xlabel('Spatial Frequency \vartheta_c');