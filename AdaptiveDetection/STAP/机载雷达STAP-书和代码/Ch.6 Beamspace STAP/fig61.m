%% Figure 61. Beamspace post-Doppler in a clutter-plus-jamming environment. K = 4. Ks = Kt = 2.
%                         Displaced Filter Beamspace post-Doppler Approach with and without 
%                         Two-Step-Nulling (TSN).
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only.
% All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters:
fo = 450*1e6;            % Operating Frequency in Hz
Pt = 200*1e3;            % Peak Transmit Power 200 kW
dc = 0.06;                   % Duty Factor or Duty Cycle 6%
Gt = 22;                       % Transmit Gain in dB
Gr = 10;                       % Column Receive Gain in dB
B  = 4*1e6;                 % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                         % System Losses in dB
fr = 300;                       % PRF in Hz
Tr = 1/fr;                       % PRI in sec.
M = 18;                        % Number of Pulses per CPI.
Tp = 200*1e-6;            % Pulse Width in sec.
%Antenna Array Parameter Value
N = 18;                         % Number of Elements
% Number of Elements z: 4
% Element Pattern: Cosine
Gel = 4;                         % Element Gain in dB
% Transmit Taper: Uniform
be = -30;                       % Element Backlobe Level in db
Nc = 360;                      % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;         % Speed of Light in m/sec.
lambda = c/fo;              % Operating wavelength in meters.
d = lambda/2;               % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:179;
Lphi = length(phi);
f = zeros(1,Lphi);
AF = zeros(1,Lphi);      % Array Factor pre-allocation.

% Platform Parameters:
beta = 1;                        % beta parameter > 1 signifies Doppler ambiguous clutter.
ha = 9000;                     % in meters.

% Noise Power Computations:
k = 1.3806488*1e-23;  % Boltzmann Constant in J/K.
To = 270;                        % Standard room Temperature in K.
F   = 3;                            % Receiver Noise Figure in dB;
Ts = To*10^(F/10);        % Receiver Temperature in K.        
Nn = k*Ts;                       % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                     % Receiver Noise Power in Watts 
sigma2 = 1;                    % Mormalized Noise Power in Watts.

% Clutter Patch Geometry computations: 
Rcik = 130000;              % (clutter) range of interest in meters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                      % Radar Range Resolution in meters.
Re = 6370000;               % Earth Radius in meters. 
ae = 4/3*Re;                   % Effective Earth Radius in meters.
% psi = -asin((Rc*Rc - ha*(ha+ 2*ae))/(2*Rc*ae));   % Grazing angle at the clutter patch in rad
                                        % (spherical earth model).
psi = asin(ha/Rcik);       % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                     % Elevation (look-down angle). Flat earth assumption.
gamma = 10^(-3/10);    % Terrain-dependent reflectivity factor. 
bwdth = 2*asin(1.4*lambda/N/d/pi);  % Antenna Azimuth 3-dB beamwidth in rad.
phia = 0;                          % Velocity Misalignment angle in degrees.

%% Clutter Covariance Matrix Computations:
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
                                                                            - sin(steering_angle*pi/180))).*cos(phi(k)*pi/180));
end

% Calculate the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Clutter Patch RCS Calculation:
PatchArea = Rcik*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

% Calculate the Clutter to Noise Ration (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);
Ksic = sigma2*diag(ksi);

% Platform Velocity for beta parameter value:
va = round(beta*d*fr/2);
Ita = d/lambda*cos(theta);   

% Calculate Spatial and Doppler Frequencies.
% Spatial frequency of the k-th clutter patch.
fspc = Ita*sin(phi*pi/180); 

% Normalized Doppler Frequency:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);

% Clutter Steering Vector Pre-allocation for sub-CPI of size K:
ac = zeros(N,Nc);
bc = zeros(M,Nc);
Vc = zeros(M*N,Nc);
       
for k=1:Nc
      ac(:,k) = exp(-1i*2*pi*fspc(k)*(0:N-1));              % Clutter Spatial Steering Vector.
      bc(:,k) = exp(-1i*2*pi*omegac(k)*(0:M-1));      % Clutter Temporal Steering Vector.   
      Vc(:,k) = kron(bc(:,k),ac(:,k));                             % Clutter Space-Time Steering Vector.
end
     
Rc = real(Vc*Ksic*Vc');     
         
Rn = sigma2*eye(M*N);
         
 %% Jammer Covariance Matrix Calculation.
     % Spatial frequency of the j-th jammer.
     J = 2;      % Number of Jammers.
     thetaj = 0; phij = [-40 25]; % Jammer elevation and azimuth angles in degrees.
     R_j = [370 370]*1e3;
     Sj = 1e-3;     % Jammer EIRP in Watts/Hz.
     fspj = d/lambda*cos(thetaj*pi/180)*sin(phij*pi/180);
     Lr = 2.0;           % System Losses on Receive in dB.
     Aj = zeros(N,J);
     for j=1:J
           Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));          % j-th Jammer Spatial Steering Vector.
     end
     
     indices = zeros(1,J);
     for j=1:J
          indices(j) = find(phi == phij(j));
     end
     grgn = grgain(indices);
     Jo = (Sj*grgn*lambda^2)./((4*pi)^2.*10^(Lr/10).*R_j.^2);
     ksi_j = Jo./Nn;
     
     Ksi_j = sigma2*diag(ksi_j);
     % Spatial Jamming Covariance Matrix
     Phi_j = Aj*Ksi_j*Aj';
     
     Rj = kron(eye(M),Phi_j);
     
%% Analytic Interference Covariance Matrix for Fully Adaptive STAP:
     Ru = Rc + Rj + Rn;    
     InvRu = inv(Ru);
     
%% Spatial jammer-plus-noise Covariance Matrix:
     Phi_jn = Phi_j + sigma2*eye(N);
     
%% Target Space-Time Steering Vector:
phit = 0; thetat = 0;                              % Target azimuth and elevation angles in degrees.
fdt = 100;                                             % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
omegat = fdt/fr;
at = exp(-1i*2*pi*fspt*(0:N-1)).';         % Target Spatial Steering Vector.
bt = exp(-1i*2*pi*fdt*(0:M-1)).';           % Target Doppler Steering Vector.
ta = chebwin(N,30);
tb = chebwin(M,30);
fd = 0:.5:300;   Lfd = length(fd);
omegad = fd/fr;
SNRo = M*N;

%% LSINR Computation for Optimum Fully Adaptive Case:
LSINRopt = zeros(1,Lfd);
for n=1:Lfd
      bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Target Doppler Steering Vector   
      vt = kron(bt,at);
      w = InvRu*vt;
      LSINRopt(n) = real(w'*vt)/SNRo;
end
    
%% Displaced Beam pre-Doppler Calculations:
Kt =  2;                                        % sub-CPI length (fixed). 
Ks = 2;                                        % sub-apperture length or number of beams used. 
M1 = M - Kt +1;                          % Number of possible sub-CPI's.
N1 = N - Ks +1;                          % Number of possible sub-appertures.
    
%% Create Doppler Filterbank Matrix F for Displaced Filter Beamspace post-Doppler STAP:
dopplerfilterbank = linspace(0,300,M+1);
omegadopplerbank = dopplerfilterbank/fr;
F0 = zeros(M1,M);
for m=1:M
      F0(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));
end
   
F30 = diag(chebwin(M1,30))*F0;                     % Doppler Filter Tapering.
        
%% Beamformer G Matrix Construction for Displaced Filter Beamspace post-Doppler STAP
%% with and without Two Step Nulling:
taper_level = 30;                                            % Spatial Tapering level in db.
g0 = exp(-1i*2*pi*fspt*(0:N1-1)).';
gtap = chebwin(N1,taper_level).*g0;           % spatial taper application.

% Calculation of Beamformer Matrix for Two-Step Nulling Approach:
GTSN = zeros(N,Ks);
for p=0:Ks-1
     Jp = [zeros(p,N1); eye(N1); zeros(N-N1-p,N1)];    % Beam Selector Matrix Jsm.
     GTSN(p+1:N1+p,p+1) = inv(Jp.'*Phi_jn*Jp)*gtap;
end

GNoTSN = toeplitz([gtap;   zeros(Ks-1,1);],  [gtap(1) zeros(1,Ks-1)]);  % N1 x Ks Beamformer matrix G.  

SINRNoTSNdf_mat  = zeros(M,Lfd);
SINRTSNdf_mat = zeros(M,Lfd);

 %% SINR Computations for Displaced Filter Beamspace post-Doppler STAP:
     % Solve a separate adaptive problem in each Doppler bin m:
for m=1:M
            
         Fm0   = toeplitz([F0(:,m);    zeros(Kt-1,1)],[F0(1,m)   zeros(1,Kt-1)]);          
         Fm30 = toeplitz([F30(:,m); zeros(Kt-1,1)],[F30(1,m) zeros(1,Kt-1)]);
         
         TmNoTSN   = kron(Fm0,GNoTSN);
         TmTSN = kron(Fm30,GTSN);
         
         Rum0   = TmNoTSN'*Ru*TmNoTSN;
         Rum30 = TmTSN'*Ru*TmTSN;
         
         bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).'; 
         gt = kron(bdfb,at);
                         
         utmNoTSM   = TmNoTSN'*gt;
         utmTSM = TmTSN'*gt;
         
         wmNoTSN   = Rum0\utmNoTSM;
         wmTSN = Rum30\utmTSM;
         
         wNoTSN   = TmNoTSN*wmNoTSN;
         wTSN = TmTSN*wmTSN;
         
        for n=1:Lfd
            bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
            vt = kron(bt,at);
            SINRNoTSNdf_mat(m,n)   = abs(wNoTSN'*vt)^2/real(wNoTSN'*Ru*wNoTSN);
            SINRTSNdf_mat(m,n) = abs(wTSN'*vt)^2/real(wTSN'*Ru*wTSN);
        end
end

SINRNoTSNdf   = max(abs(SINRNoTSNdf_mat));
SINRTSNdf = max(abs(SINRTSNdf_mat));
   
LSINRNoTSNdf     = SINRNoTSNdf/SNRo;
LSINRTSNdf   = SINRTSNdf/SNRo;
     
  %% Plot the SINR Losses:
     figure('NumberTitle', 'off','Name', 'Figure 61. Beamspace post-Doppler in a clutter-plus-jamming environment. K=4, Ksm = Ktm = 2.');
    
 %  Displaced Filter Beamspace post-Doppler STAP Results:
     plot(fd,10*log10(LSINRopt),'k','LineWidth',1.5)      
     hold on;
     plot(fd,10*log10(LSINRNoTSNdf),'LineWidth',1.5)      
     plot(fd,10*log10(LSINRTSNdf),'g','LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
      title('Ks = Kt = 2');
     ylim([-30 1]);
     xlim([-5 305]);
     legend('Optimum','No TSN','With TSN','Location','South');      
    