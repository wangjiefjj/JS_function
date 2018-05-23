%% Figure 59. Beamspace post-Doppler in a clutter-only environment. K = 9. Ks = Kt = 3.
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
         
%% Analytic Interference Covariance Matrix for Fully Adaptive STAP:
Ru = Rc + Rn;    
InvRu = inv(Ru);
     
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
Kt =  3;                                        % sub-CPI length (fixed). 
Ks = 3;                                        % sub-apperture length or number of beams used. 
M1 = M - Kt +1;                          % Number of possible sub-CPI's.
N1 = N - Ks +1;                          % Number of possible sub-appertures.
    
%% Create Doppler Filterbank Matrix F for Displaced Filter Beamspace post-Doppler STAP:
dopplerfilterbank = linspace(0,300,M+1);
omegadopplerbank = dopplerfilterbank/fr;
F0 = zeros(M1,M);
for m=1:M
      F0(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));
end
   
F40 = diag(chebwin(M1,40))*F0;
        
%% Beamformer G Matrix Construction for Displaced Filter Beamspace post-Doppler STAP:
g0 = exp(-1i*2*pi*fspt*(0:N1-1)).';
g30 = chebwin(N1,30).*g0;           % 30-dB spatial taper.
G0   = toeplitz([g0;   zeros(Ks-1,1);],  [g0(1) zeros(1,Ks-1)]);  % N1 x Ks Beamformer matrix G.  
G30 = toeplitz([g30; zeros(Ks-1,1);], [g30(1) zeros(1,Ks-1)]);  % N1 x Ks Beamformer matrix G.  

SINR0df_mat  = zeros(M,Lfd);
SINR30df_mat = zeros(M,Lfd);
SINR0af_mat  = zeros(M,Lfd);
SINR30af_mat = zeros(M,Lfd);

%%  Create Doppler Filterbank Matrix F for Adjacent Filter Beamspace post-Doppler STAP:
U2 = zeros(M,M);
Pt = floor(Kt/2);
if  mod(Kt,2)    % If Kt is odd
    for m=1:M
         U2(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1));     % Doppler Filter Steering Vector   
    end
else            % if Kt is even:
      outomegadoppler = zeros(1,M);
     for m=1:M
          outomegadoppler(m) = (omegadopplerbank(m) + omegadopplerbank(m+1))/2;
          U2(:,m) = 1/sqrt(M)*exp(-1i*2*pi*outomegadoppler(m)*(0:M-1));     % Doppler Filter Steering Vector   
    end
end
       
td0af   = ones(M,1);                                                       % Uniform Doppler Taper.
td40af = chebwin(M,40);                                               % 40-dB Chebyshev Doppler Taper.
     
% Create Doppler Filter Bank in Fab matrix for Adjacent Filter Beamspace post-Doppler method:
Faf0   = diag(td0af)*U2;                               
Faf40 = diag(td40af)*U2;
             
Fmaf0   = zeros(M,Kt,M);
Fmaf40 = zeros(M,Kt,M);
  
for m=1:M
     if mod(Kt,2)                  % if K is odd and >1.
        if (m-Pt>0) && (m+Pt<=M)
            Fmaf0(:,:,m)   = Faf0(:,m-Pt:m+Pt);                                                          % Eq. 231.
            Fmaf40(:,:,m) = Faf40(:,m-Pt:m+Pt);
        elseif (m-Pt<=0) && (m+Pt<=M)
              Fmaf0(:,:,m)   = [Faf0(:,M+(m-Pt):M)      Faf0(:,1:m+Pt)];                      % Eq. 231.
              Fmaf40(:,:,m) = [Faf40(:,M+(m-Pt):M) Faf40(:,1:m+Pt)];       
        elseif m+Pt>M
              Fmaf0(:,:,m)   = [Faf0(:,m-Pt:M)      Faf0(:,1:m+Pt-M)];                          % Eq. 231.
              Fmaf40(:,:,m) = [Faf40(:,m-Pt:M) Faf40(:,1:m+Pt-M)];                                       
        end
                      
     else      % if K is even.
              
        if (m-Pt>0) && (m+Pt<=M+1)
            Fmaf0(:,:,m)   = Faf0(:,m-Pt:m+Pt-1);                                                % Eq. 231.
            Fmaf40(:,:,m) = Faf40(:,m-Pt:m+Pt-1);
        elseif (m-Pt<=0) && (m+Pt<=M)
            Fmaf0(:,:,m)   = [Faf0(:,M+(m-Pt):M)      Faf0(:,1:m+Pt-1)];                % Eq. 231.
            Fmaf40(:,:,m) = [Faf40(:,M+(m-Pt):M) Faf40(:,1:m+Pt-1)];       
        elseif m+Pt>M+1
             Fmaf0(:,:,m)   = [Faf0(:,m-Pt:M)        Faf0(:,1:m-M+Pt-1)];           % Eq. 231.
             Fmaf40(:,:,m) = [Faf40(:,m-Pt:M)   Faf40(:,1:m-M+Pt-1)];                                       
        end 
    end
                  
end

%%  Create Beamformer Matrix G for Adjacent Filter Beamspace post-Doppler STAP:
selbeamdist = 6.35;                     % Selected Beam Angular Distance in degrees.
beamangles1 = [-9:-1  0:8]*selbeamdist; thetabeam = theta;
  
beamangleinc = beamangles1(2) - beamangles1(1);
beamangles2 = beamangles1 + beamangleinc/2;
     
beamfreqs1 = d/lambda*cos(thetabeam*pi/180)*sin(beamangles1*pi/180);
beamfreqs2 = d/lambda*cos(thetabeam*pi/180)*sin(beamangles2*pi/180);
     
Godd   = zeros(N,N);
Geven = zeros(N,N);
     
for n=1:N
     Godd(:,n)  = 1/sqrt(N)*exp(-1i*2*pi*beamfreqs1(n)*(0:N-1));
     Geven(:,n) = 1/sqrt(N)*exp(-1i*2*pi*beamfreqs2(n)*(0:N-1));
end

%% Create Beam Selector Matrix Jsm for Adjacent Filter Beamspace post-Doppler STAP:
Ps = floor(Ks/2);
if mod(Ks,2)  % if Ks is odd
     Jsm = [zeros(N/2-Ps,Ks); eye(Ks); zeros(N/2-Ps-1,Ks)];    % Beam Selector Matrix Jsm.
    G0af = Godd;                                                                          % N x Ks Beamformer matrix.
  G30af = diag(chebwin(N,30))*G0af;
else
     Jsm = [zeros(N/2-Ps,Ks); eye(Ks); zeros(N/2-Ps,Ks)];        % Beam Selector Matrix Jsm.    
    G0af = Geven;                                                                          % N x Ks Beamformer matrix.
  G30af = diag(chebwin(N,30))*G0af;
end

G0maf = G0af*Jsm;
G30maf = G30af*Jsm;

     % Solve a separate adaptive problem in each Doppler bin m:
     for m=1:M
         
         %% A. SINR Computations for Displaced Filter Beamspace post-Doppler STAP:
        
         Fm0   = toeplitz([F0(:,m);    zeros(Kt-1,1)],[F0(1,m)   zeros(1,Kt-1)]);          
         Fm30 = toeplitz([F40(:,m); zeros(Kt-1,1)],[F40(1,m) zeros(1,Kt-1)]);
         
         Tm0   = kron(Fm0,G0);
         Tm30 = kron(Fm30,G30);
         
         Rum0   = Tm0'*Ru*Tm0;
         Rum30 = Tm30'*Ru*Tm30;
         
         bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).'; 
         gt = kron(bdfb,at);
                         
         utm0   = Tm0'*gt;
         utm30 = Tm30'*gt;
         
         wm0   = Rum0\utm0;
         wm30 = Rum30\utm30;
         
         w0   = Tm0*wm0;
         w30 = Tm30*wm30;
         
        for n=1:Lfd
            bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
            vt = kron(bt,at);
            SINR0df_mat(m,n)   = abs(w0'*vt)^2/real(w0'*Ru*w0);
            SINR30df_mat(m,n) = abs(w30'*vt)^2/real(w30'*Ru*w30);
        end
     end

SINR0df   = max(abs(SINR0df_mat));
SINR30df = max(abs(SINR30df_mat));
   
LSINR0df     = SINR0df/SNRo;
LSINR30df   = SINR30df/SNRo;
     
%% B. SINR Computations for Adjacent Filter Beamspace post-Doppler STAP:
for m=1:M
   
     F0maf   = Fmaf0(:,:,m);                                   
     F30maf = Fmaf40(:,:,m);
       
     Tm0af = kron(F0maf,G0maf);
     Tm30af = kron(F30maf,G30maf);
       
      Rum0af   = Tm0af'*Ru*Tm0af;
      Rum30af = Tm30af'*Ru*Tm30af;
         
      bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).'; 
      gt = kron(bdfb,at);
         
      utm0af   = Tm0af'*gt;
      utm30af = Tm30af'*gt;
         
      wm0af   = Rum0af\utm0af;
      wm30af = Rum30af\utm30af;
         
      w0af   = Tm0af*wm0af;
      w30af = Tm30af*wm30af;
              
      for n=1:Lfd
            bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
            vt = kron(bt,at);
            SINR0af_mat(m,n)   = abs(w0af'*vt)^2/real(w0af'*Ru*w0af);
            SINR30af_mat(m,n) = abs(w30af'*vt)^2/real(w30af'*Ru*w30af);
      end
             
end
   
SINR0af   = max(abs(SINR0af_mat));
SINR30af = max(abs(SINR30af_mat));
   
LSINR0af     = SINR0af/SNRo;
LSINR30af   = SINR30af/SNRo;

  %% Plot the SINR Losses:
     figure('NumberTitle', 'off','Name', 'Figure 59. Beamspace post-Doppler in a clutter only environment. K=9, Ksm = Ktm = 3.');
    
 %  Displaced Filter Beamspace post-Doppler STAP Results:
     plot(fd,10*log10(LSINRopt),'k','LineWidth',1.5)      
     hold on;
     plot(fd,10*log10(LSINR0df),'LineWidth',1.5)      
     plot(fd,10*log10(LSINR30df),'g','LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 1]);
     xlim([-5 305]);
          
 %  Adjacent Filter Beamspace post-Doppler STAP Results:
     plot(fd,10*log10(LSINR0af),'r','LineWidth',1.5)      
     hold on;
     plot(fd,10*log10(LSINR30af),'m','LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 1]);
     xlim([-5 305]);
     legend('Optimum','Displaced-Uniform', 'Displaced-Tapered', 'Adjacent-Uniform', 'Adjacent-Tapered','Location','South');      