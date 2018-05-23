%% Figure 14. Example clutter ridges with velocity misalignment, for Doppler-ambiguous clutter.

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
d = lambda/2;                 % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:0.5:180;
Lphi = length(phi);
fd = zeros(1,Lphi);

%% Platform Parameters.
beta = 1.5:0.5:3;             % Beta parameter values.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Depression angle to ik-th clutter patch (flat earth model).
phia = [30 45 60 75];         % Velocity Misalignment angle in degrees.

%% Clutter Spatial Frequency of ik-th Clutter Patch.
fsp = d/lambda*cos(theta)*sin(phi*pi/180);

%% Plot Normalized Doppler Frequency vs Spatial Frequency.
figure('NumberTitle', 'off','Name', ...
    'Figure 14. Example clutter ridges with velocity misalignment, for Doppler-ambiguous clutter.', ...
    'Position', [1 1 1150 1250]);

for i1=1:length(beta)
    for i2 = 1:length(phia)
        % Platform Velocity for various beta parameter values:
        va = beta(i1)*d*fr/2;
        
        % Doppler Frequency from ik-th Clutter Patch
        fd = 2*va/lambda*cos(theta)*sin(phi*pi/180 + phia(i2)*pi/180);      % Eq. (82)
        
        % Normalized Doppler Frequency:
        omegac = fd*Tr;
        
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
        
        i3 = 4*(i1-1)+i2;
        subplot(4,4,i3);
        plot(fsp,front,'.');
        hold on;
        plot(fsp,back,'r.');
        title(['\beta = ', num2str(beta(i1)), ',  \phi_a = ',num2str(phia(i2)),'\circ']);
        if i3==1||i3==5||i3==9||i3==13
            ylabel('Norm. Doppler Freq. \omega_c');
        end
        ylim([-0.5 0.5]); xlim([-0.5 0.5]);
        grid on;
        if i3==16
            legend('front','back','Location','Best');
        end
        if i1==4
            xlabel('Spatial Frequency \vartheta_c');
        end
    end
end

tightfig;
