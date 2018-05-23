%% Figure 29. Expected SINR Loss for SMI with matched steering vector.
%%
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only. All rights reserved.

clc;  clear; close all;

%% Radar System Operational Parameters
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
M = 18;                       % Number of Pulses per CPI.
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

% Platform Parameters:
beta = 1;                     % beta parameter.
ha = 9e3;                     % Platform altitude in meters.

%% Thermal Noise Power Computations
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power in Watts.

%% Clutter Patch Geometry computations
Rcik = 130000;                % (clutter) range of interest in meters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Elevation (look-down angle) in rad. Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
phia = 0;                     % Velocity Misalignment angle in degrees.

%% Clutter-to-Noise Ratio (CNR) Calculation
% Calculate the Voltage Element Pattern:
for i =1:Lphi
    if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
    else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
    end
end

% Calculate the Array Factor (AF) (Voltage):
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
    AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ...
        - sin(steering_angle*pi/180))));
end

% Calculate the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Clutter Patch RCS Calculation:
PatchArea = Rcik*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

% Calculate the Clutter to Noise Ratio (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);
Ksic = sigma2*diag(ksi);

%% Clutter Covariance Matrix Computations
% Platform Velocity for beta parameter value:
va = round(beta*d*fr/2);
Ita = d/lambda*cos(theta);

% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
% Spatial frequency of the k-th clutter patch:
fsp = Ita*sin(phi*pi/180);
% Normalized Doppler Frequency of the k-th clutter patch:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);

% Clutter Steering Vector Pre-allocation:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1));    % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));           % Space-Time Steering Vector.
end

Rc = Vc*Ksic*Vc';                            % Eq. (64)

Rn = sigma2*eye(M*N);

%% Jamming Covariance Matrix Calculation
J = 2;                                                    % Number of Jammers.
thetaj = 0; phij = [-40 25];                              % Jammer elevation and azimuth angles in degrees.
R_j = [370 370]*1e3;
Sj = 1e-3;                                                % Jammer ERPD in Watts/Hz.
fspj = d/lambda*cos(thetaj*pi/180)*sin(phij*pi/180);      % Spatial frequency of the j-th jammer.
Lrj = 1.92;                                               % System Losses on Receive in dB.
Aj = zeros(N,J);
for j=1:J
    Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));              % Jammer Spatial Steering Vector.
end

indices= zeros(1,J);
for j=1:J
    indices(j) = find(phi == phij(j));
end
grgn = grgain(indices);
ksi_j = (Sj*grgn*lambda^2)./((4*pi)^2.*Nn*10^(Lrj/10).*R_j.^2);

Ksi_j = sigma2*diag(ksi_j);
Phi_j = Aj*Ksi_j*Aj';                                     % Eq. (47)
% Jamming Covariance Matrix:
Rj = kron(eye(M),Phi_j);                                  % Eq. (45)

%% Analytic Interference Covariance Matrix Calculation
Ru = Rc + Rj + Rn;

%% SINR Loss calculations:
fd = -150:.5:150;   Lfd = length(fd);
phit = 0; thetat = 0;                                     % Target Azimuth and Elevation Angles.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);      % Target Spatial Frequency.
omegat = fd/fr;
at = exp(1i*2*pi*fspt*(0:N-1));                           % Target Spatial Steering Vector.
LSINRopt = zeros(1,Lfd);
InvRu = inv(Ru);
SNRo = M*N;                                               % Eq. (119)
Ndof = M*N;                                               % Number of adjustable weights (Degrees of Freedom)
Samples = [Ndof 2*Ndof 3*Ndof 4*Ndof 5*Ndof];             % Number of Samples(n3) used for SMI.
SINRopt = zeros(1,Lfd);

% Calculate LSINRopt:
for n1=1:Lfd
    bt = exp(1i*2*pi*omegat(n1)*(0:M-1));                 % Dummy Target Doppler Steering Vector
    vt = kron(bt,at).';
    wopt = InvRu*vt; %#ok<*MINV>
    SINRopt(n1) = real(wopt'*vt);                         % Eq. (112)
    LSINRopt(n1) = SINRopt(n1)/SNRo;                      % Eq. (120) 
end

%% Calculate by Monte Carlo Experiments and Plot the Expected SINR Loss 
figure('NumberTitle', 'off','Name', 'Figure 29. SINR Loss for SMI with matched steering vector',...
       'Position', [1 1 900 700]);
plot(fd,10*log10(abs(LSINRopt)),'LineWidth',1.5)
hold on;

% First Calculate LSINRest expected value using Monte Carlo method:
NRuns = 40;                         % Number of Monte Carlo runs.
rho = zeros(NRuns,Lfd);
SINRa = zeros(NRuns,Lfd);
time_loop = zeros(1,NRuns);
% Set the seed of the RNG.
LSINRest = zeros(Lfd,length(Samples));
colors = [0 1 0; 1 0 0 ; 1 1 0; 0 1 1; 1 0 1;];
rng(131);
warning('off', 'all');

for n3=1:length(Samples)
    
    for n2=1:NRuns
        % n2 %#ok<*NOPTS>
        X   = zeros(M*N,Samples(n3));
        Rest = zeros(M*N,M*N);
        for n=1:Samples(n3)
            
            % Create Thermal Noise Measurement space-time vector:
            chi_n = sqrt(sigma2/2)*(randn(M*N,1) + 1i*randn(M*N,1));
            
            % Create #1 Jammer Interference Measurement temporal vector:
            alphaj1 = sqrt(sigma2*ksi_j(1,1)/2)*(randn(M,1) + 1i*randn(M,1));
            % Create #2 Jammer Interference Measurement temporal vector:
            alphaj2 = sqrt(sigma2*ksi_j(1,2)/2)*(randn(M,1) + 1i*randn(M,1));
            % Create #1 Jammer Interference Measurement space-time vector:
            chi_j1 = kron(alphaj1, Aj(:,1));
            % Create #2 Jammer Interference Measurement space-time vector:
            chi_j2 = kron(alphaj2, Aj(:,2));
            
            % Create Total Jamming Interference Measurement space-time vector:
            chi_j = chi_j1 + chi_j2;
            
            % Create Total Clutter Interference Measurement space-time vector:
            R = rand(1,Nc);  I = randn(1,Nc);
            Ksi1 = repmat(sqrt(sigma2*ksi/2).*(R + 1i*I),[M*N 1]);
            chi_ik = sum(Ksi1.*Vc,2);
            % Add the interference component vectors to get the total
            % (target free) measured interference space-time snapshot:
            X(:,n) = chi_ik + chi_j + chi_n;
            
        end
        
        % Interference Sample Covariance Matrix Computation:
        meanX = 1/Samples(n3)*sum(X,2);
        Ruest = 1/Samples(n3)*(X*X') - (meanX*meanX');                         % Eq. (129)
        InvRuest = inv(Ruest);
        
        for n1=1:Lfd
            bt = exp(1i*2*pi*omegat(n1)*(0:M-1)); % Dummy Target Doppler Steering Vector
            vt = kron(bt,at).';
            w = InvRu*vt;
            west = InvRuest*vt;                                                % Eq. (130)
            rho(n2,n1) = abs(west'*vt)^2/real(west'*Ru*west)/real(w'*vt);      % Eq. (132)
        end
    end
    
    % Average the Estimated SINR Loss Factor ñ (rho) and apply it to the Optimal SINR Loss
    % LSINTopt to calculate the estimated SMI SINR Loss:
    LSINRest(:,n3) = 1/NRuns*sum(rho).*LSINRopt;                               % Eq. (136)
    
    plot(fd,10*log10(abs(LSINRest(:,n3))),'LineWidth',1.5,'Color', colors(n3,:))
    
end

ylabel('Expected SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
legend('Known Covariance Matrix', 'SMI using Ke = N_{dof} samples', ...
    'SMI using Ke = 2N_{dof} samples', 'SMI using Ke = 3N_{dof} samples', ...
    'SMI using Ke = 4N_{dof} samples','SMI using Ke = 5N_{dof} samples', ...
    'Location','East')
grid on;
