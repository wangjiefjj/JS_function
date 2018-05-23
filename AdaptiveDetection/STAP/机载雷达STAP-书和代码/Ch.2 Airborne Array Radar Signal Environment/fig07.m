%% Figure 7. Clutter Loci for Different Platform Velocities.
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
c = 299792458;                % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing.

% Azimuth angle in degrees:
phi = -180:.5:180;
Lphi = length(phi);

%% Platform Parameters.
va = 50;                      % Platform velocity in m/sec.
ha = 9e3;                     % Platform altitude in meters.
beta = 0:0.5:7;               % Beta parameter values.
Rc = 13e4;                    % (clutter) range of interest in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Depression angle to ik-th clutter patch (flat earth model).

% Clutter Spatial Frequency of ik-th Clutter Patch:
fsp = d/lambda*cos(theta)*sin(phi*pi/180);

%% Calculate and Plot the Clutter Loci for Different Platform Velocities.
figure('NumberTitle', 'off','Name', ...
       'Figure 7. Clutter Loci for Different Velocities of Side Looking Airborne Radar (SLAR) ', ...
       'Position', [50 50 1150 700] );

for i=1:length(beta)
    % Platform Velocity for various beta parameter values:
    va = beta(i)*d*fr/2;                                   % Eq. (71)
    
    % Doppler Frequency from ik-th Clutter Patch
    % fd = 2*va/lambda*cos(theta)*sin(phi*pi/180);
    
    % Normalized Doppler Frequency:
    % omegac = 2*va*Tr/d*fsp;
    omegac = beta(i)*fsp;                                  % Eq. (70)
    
    % This loop simulates the fold-over of the Clutter ridge (i.e. brings omegac 
    % back into the [-0.5 0.5] interval) when clutter is Doppler-ambiguous (â>1).
    for k=1:Lphi
        if omegac(k)>0.5 && omegac(k)<=1.5
            omegac(k) = omegac(k) - 1;
        elseif omegac(k)<-0.5 && omegac(k)>=-1.5
            omegac(k) = omegac(k) + 1;
        end
        
        if omegac(k)>1.5 && omegac(k)<=2.5
            omegac(k) = omegac(k) - 2;
        elseif omegac(k)<-1.5 && omegac(k)>=-2.5
            omegac(k) = omegac(k) + 2;
        end
        
        if omegac(k)>2.5 && omegac(k)<=3.5
            omegac(k) = omegac(k) - 3;
        elseif omegac(k)<-2.5 && omegac(k)>=-3.5
            omegac(k) = omegac(k) + 3;
        end
    end
    
    % Plot Normalized Doppler Frequency vs Spatial Frequency.
    subplot(3,5,i);
    plot(fsp,omegac,'.');
    title(['v_a = ',num2str(round(va)),' m/sec, \bf\beta = ',num2str(beta(i))]);
    if i==1 || i== 6 || i==11
        ylabel('Norm. Doppler Frequency \omega_c');
    end
    if i>10
        xlabel('Spatial Frequency \vartheta_c');
    end
    ylim([-0.5 0.5]); xlim([-0.5 0.5]);
    grid on;
    
end

tightfig;
