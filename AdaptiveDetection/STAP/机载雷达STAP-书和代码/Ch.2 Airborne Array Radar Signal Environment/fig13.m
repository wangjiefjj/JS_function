%% Figure 13. Example Clutter Ridges with Velocity Misalignment, for â=1.
%%
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters.
fo = 450e6;                   % Operating Frequency in Hz
fr = 300;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
c = 299792458;                % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing in meters.

% Azimuth angle in degrees:
phi = -180:0.25:180;
Lphi = length(phi);
fd = zeros(1,Lphi);

%% Platform Parameters.
beta = 1;                     % Beta parameter value.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.
va = 50;                      % Platform Velocity on m/sec.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Depression angle to ik-th clutter patch (flat earth model).
phia = 0:15:105;              % Velocity Misalignment angle in degrees.

% Clutter Spatial Frequency of ik-th Clutter Patch:
fsp = d/lambda*cos(theta)*sin(phi*pi/180);

%% Plot Normalized Doppler Frequency vs Spatial Frequency.
figure('NumberTitle', 'off','Name', ...
    'Figure 13. Example Clutter Ridges with Velocity Misalignment, for â=1', 'Position', [50 50 1150 550]);

for i=1:length(phia)
    % Platform Velocity for various beta parameter values:
    % Doppler Frequency from ik-th Clutter Patch
    fd = 2*va/lambda*cos(theta)*sin(phi*pi/180 + phia(i)*pi/180);   % Eq. (82)
    
    % Normalized Doppler Frequency:
    omegac = fd*Tr;
    
    front = zeros(1,Lphi);
    back = zeros(1,Lphi);
    
    for k=1:Lphi
        if abs(phi(k)) <= 90
            front(k) = omegac(k);
        else
            back(k) = omegac(k);
        end
    end
    
    zeroindfront = find(front==0);
    zeroindback = find(back==0);
    front(zeroindfront)   = NaN;
    back(zeroindback) = NaN;
    
    subplot(2,4,i);
    plot(fsp,front,'.');
    hold on;
    if i>1
        plot(fsp,back,'r.');
    end
    
    title(['\phi_a = ',num2str(phia(i)),'\circ']);
    if i==1 || i==5
       ylabel('Norm. Doppler Frequency \omega_c');
    end
    if i>4
      xlabel('Spatial Frequency \vartheta_c');
    end
    ylim([-0.5 0.5]); xlim([-0.5 0.5]);
    grid on;
    if i==8
        legend('front','back','Location','Best');
    end
end

tightfig;
