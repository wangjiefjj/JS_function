%% Radar System Parameter
fo = 450e6;                 %Operating Frequency
Pt = 200e3;                 %Peak Power
Df = 0.06;                  %Duty factor
Gt = 22;                    %Transmit Gain dB
Gr = 10;                    %Column Receive Gain dB
B = 4e6;                    %Instantaneous Bandwidth
Nf = 3;                     %Noise Figure dB
Ls = 4;                     %System Losses dB
fr = 300;
Tr = 1/fr;
M = 18;                     %Number of Pulses/CPI
c = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
Tp = 200e-6;                %Pulse Width
% Azimuth angle in degrees:
phi = -180:.5:180;
Lphi = length(phi);
%% Anrenna Array Parameter
N = 18;                     %Number of Elements
Gel = 4;                    %Element Gain dB
d = lambda/2;                 % Interelement Spacing.
%% Platform Parameters:
va = 50;                      % Platform velocity in m/sec.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Depression angle to ik-th clutter patch (flat earth model).

