%% Figure 68.(a) SINR loss for example system with -40dB backlobe, â=2.6, ó_õ = 0.2 m/sec, 
%%                                                                                                                 90 deg vel. misalinment.
%
%
% Coded by Ilias Konsoulas, 25 Mar. 2016.
% Code provided for educational purposes only.
% All rights reserved.

clc;  clear; close all;

%% Radar System Operational Parameters:
fo = 450e6;              % Operating Frequency in Hz
Pt = 200e3;              % Peak Transmit Power 200 kW
dc = 0.06;                  % Duty Factor or Duty Cycle 6%
Gt = 22;                     % Transmit Gain in dB
Gr = 10;                     % Column Receive Gain in dB
B  = 4e6;                   % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                        % System Losses in dB
fr = 300;                      % PRF in Hz
Tr = 1/fr;                      % PRI in sec.
M = 18;                       % Number of Pulses per CPI.
K = 3;                          % Number of Pulses per sub-CPI.
M2 = M - K + 1;          % Number of sub-CPI's for K=3.
Tp = 200e-6;              % Pulse Width in sec.

%Antenna Array Parameter Value
N = 18;                         % Number of Elements
% Number of Elements z: 4
% Element Pattern: Cosine
Gel = 4;                         % Element Gain in dB
% Transmit Taper: Uniform
be = -40;                       % Element Backlobe Level in db
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
beta = 2.6;                     % beta parameter.
ha = 9000;                     % in meters.

% Noise Power Computations:
k = 1.3806488e-23;      % Boltzmann Constant in J/K.
To = 270;                        % Standard room Temperature in K.
F   = 3;                            % Receiver Noise Figure in dB;
Ts = To*10^(F/10);        % Receiver Temperature in K.        
Nn = k*Ts;                       % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                     % Receiver Noise Power in Watts 
sigma2 = 1;                    % Mormalized Noise Power in Watts.

% Clutter Patch Geometry computations: 
Rcik = 13e4;                   % (clutter) range of interest in meters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                      % Radar Range Resolution in meters.
Re = 637e4;                   % Earth Radius in meters. 
ae = 4/3*Re;                   % Effective Earth Radius in meters.
% psi = -asin((Rc*Rc - ha*(ha+ 2*ae))/(2*Rc*ae));   % Grazing angle at the clutter patch in rad
                                         % (spherical earth model).
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                      % Elevation (look-down angle). Flat earth assumption.
gamma = 10^(-3/10);     % Terrain-dependent reflectivity factor. 
bwdth = 2*asin(0.446*lambda/N/d);  % Antenna Azimuth 3-dB beamwidth in rad.
phia = 90;                         % Velocity Misalignment angle in degrees.
sigma_icm = 0.2;            % Intrinsic Clutter Motion Standard Deviation in meters/sec. 
kc = 4*pi*sigma_icm/lambda;

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
fsp = Ita*sin(phi*pi/180); 
% Normalized Doppler Frequency:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);
     
% Clutter Steering Vector Pre-allocation for sub-CPI of size K:
a = zeros(N,Nc);
b = zeros(M,Nc);
b2 = zeros(K,Nc);
Rc = zeros(M*N,M*N);
Rcsub2 = zeros(K*N,K*N);

gammac         =  exp(-(kc^2*Tr^2*(0:M-1).^2)/2);
gammac_sub =  exp(-(kc^2*Tr^2*(0:K-1).^2)/2);
Gammac         = toeplitz(gammac);
Gammac_sub = toeplitz(gammac_sub);

for k=1:Nc
     a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1)); % Spatial Steering Vector.
     b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Time Steering Vector   
     b2(:,k) = exp(1i*2*pi*omegac(k)*(0:K-1));  % Time Steering Vector for K = 3.   
     Rc         = Rc          + ksi(k)*kron(Gammac.*(b(:,k)*b(:,k)'), a(:,k)*a(:,k)');       
     Rcsub2 = Rcsub2 + ksi(k)*kron(Gammac_sub.*(b2(:,k)*b2(:,k)'), a(:,k)*a(:,k)');       
end

Rn = sigma2*eye(M*N);
Rnsub2 = sigma2*eye(K*N);
     
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
           Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));          % Jammer Spatial Steering Vector.
     end
     
     indices = zeros(1,J);
     for j=1:J
          indices(j) = find(phi == phij(j));
     end
     grgn = grgain(indices);
     Jo = (Sj*grgn*lambda^2)./((4*pi)^2.*10^(Lr/10).*R_j.^2);
     ksi_j = Jo./Nn;
     
     Ksi_j = sigma2*diag(ksi_j);
     Phi_j = Aj*Ksi_j*Aj';
     
     Rj = kron(eye(M),Phi_j);
     
     Rjsub2 = kron(eye(K),Phi_j);
     
     %% Spatial jammer-plus-noise Covariance Matrix:
      Phi_jn = Phi_j + sigma2*eye(N);
  
     %% Analytic Interference Covariance Matrix for Fully Adaptive STAP:
     Ru = Rc + Rj + Rn;
     
     %% 1. SINR Computations for Optimum Fully Adaptive STAP:
     fspt = 0;
     at = exp(1i*2*pi*fspt*(0:N-1)).';           % Dummy Target Spatial Steering Vector.

     fd = 0:.5:300;   Lfd = length(fd);
     omegad = fd/fr;
     dopplerfilterbank2 = linspace(0,300,M2+1);
     omegadopplerbank2 = dopplerfilterbank2/fr;
     LSINRopt = zeros(1,Lfd);
     SINRsub2_mat = zeros(length(dopplerfilterbank2),Lfd);     
     InvRu = inv(Ru);
    SNRo = M*N;
     
     % LSINR Computation for Optimum Fully Adaptive Case:
     for n=1:Lfd
          bt = exp(1i*2*pi*omegad(n)*(0:M-1)).'; % Target Doppler Steering Vector   
          vt = kron(bt,at);
          w = InvRu*vt; %#ok<MINV>
          LSINRopt(n) = real(w'*vt)/SNRo;
     end
       
     
 %% Analytic Interference Covariance Matrix for sub-CPI (K=3):
     Rusub2 = Rcsub2 + Rjsub2 + Rnsub2;
     InvRusub2 = inv(Rusub2);
      
   % Target Space-Time Steering Vector:
     phit = 0; thetat = 0;                             % Target azimuth and elevation angles in degrees.
     fdt = 100;                                             % Target Doppler Frequency.
     fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
     omegact = fdt/fr;
     at = exp(1i*2*pi*fspt*(0:N-1)).';           % Target Spatial Steering Vector.
     ta = chebwin(N,30);                              % 30 dB Chebychev Spatial Tapper. 
     
     bt2 = [1; -1; 1;];
     tb2 = [1;  2; 1;];                                       % Binomial Doppler Taper for K = 3.
     
     gt2 = kron(tb2.*bt2,ta.*at);    
     wsub2 = Rusub2\gt2;            % Weight Vector for K = 3;
    
     W2 = zeros(M*N,M2);
  % W2 matrix construction. Equation (180).
  % Assume that weight vectors for each sub-CPI are equal.
     for k=1:M2
          W2((k-1)*N+1:(k+2)*N,k) = wsub2;
     end
         
     %% 2. SINR Computations for Element-Space pre-Doppler, sub-CPI with K=3. 
      td2 =  chebwin(M2,30);                       % 30 dB Chebychev Doppler Taper.  
      for m=1:length(dopplerfilterbank2)
                um = exp(1i*2*pi*omegadopplerbank2(m)*(0:M2-1)).';   % Doppler Filter Steering Vector   
                fm2 = td2.*um;
                wm2 = W2*fm2;                                        % m-th Doppler bin composite weight vector.
                                   
           for n=1:Lfd
                 bt = exp(1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
                 vt = kron(bt,at);
                SINRsub2_mat(m,n) = abs(wm2'*vt)^2/real(wm2'*Ru*wm2);
           end
     end
     
     SINRsub2 = max(abs(SINRsub2_mat));
     LSINRsub2 = SINRsub2/SNRo;
     
     
   %% 3. Element Space post-Doppler: PRI-Staggered Post-Doppler method
    % Doppler Filter Bank Creation:
     dopplerfilterbank = linspace(0,300,M+1);
     omegadopplerbank = dopplerfilterbank/fr;
     fd = 0:.5:300;   Lfd = length(fd);
     omegad = fd/fr;
     SNRo = M*N;
     
     K = 3;
     P = floor(K/2);
     M1= M - K +1;
     
     U1 = zeros(M1,M);
     for m=1:M
          U1(:,m) =  1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));     % Doppler Filter Steering Vector   
     end

     td30 = chebwin(M1,30);                                               % 30-dB Chebyshev Doppler Taper.
     F30 = diag(td30)*U1;
     
 % Create Doppler Filter Bank in Fm Matrix for PRI-Staggered Post-Doppler method:
    Fm30 = zeros(M,K,M);
    for m=1:M
          Fm30(:,:,m) = toeplitz([F30(:,m); zeros(K-1,1)],[F30(1,m) zeros(1,K-1)]);
    end
     
    SINR30_mat     = zeros(M,Lfd);
    SINRab30_mat = zeros(M,Lfd);
    
  % PRI-Staggered SINR Computations:   
   for m=1:M
              
        f30m = Fm30(:,:,m);
        R30um = kron(f30m,eye(N))'*Ru*kron(f30m,eye(N));
        bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).'; 
        gt = kron(bdfb,ta.*at);                   % Desired Vector common for both methods.
        gt30m =  kron(f30m,eye(N))'*gt;
        w30m = R30um\gt30m;                % Calculate K*N X 1 Adaptive Weight for m-th Doppler Bin.
        w30 = kron(f30m,eye(N))*w30m;
                             
        for n=1:Lfd
              bt = exp(1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
              vt = kron(bt,at);
              SINR30_mat(m,n) = abs(w30'*vt)^2/real(w30'*Ru*w30);
        end
     
   end
     
  % PRI-Staggered SINR Loss Computations:
    SINR30 = max(abs(SINR30_mat));
    LSINR30 = SINR30/SNRo;
    
     
    %% 4. Beamspace pre-Doppler: Displaced Beams with Two Step Nulling:
    Kt =  3;                                        % sub-CPI length (fixed). 
    Ks = 5;                                        % various sub-apperture lengths or number of beams used. 
    M1 = M - Kt +1;                          % Number of possible sub-CPI's.

    % Create Doppler Filterbank Matrix F:
     dopplerfilterbank = linspace(0,300,M1+1);
     omegadopplerbank = dopplerfilterbank/fr;
     F0 = zeros(M1,M1);
     for m1=1:M1
          F0(:,m1) = 1/sqrt(M1)*exp(-1i*2*pi*omegadopplerbank(m1)*(0:M1-1));
     end
   
    F30 = diag(chebwin(M1,30))*F0;
        
    % Beamformer G Matrix Construction for Displaced Beam pre-Doppler.
     tb = [1;  2; 1;];                                         % Temporal Binomial Taper.
     bt = [1; -1; 1;];                                         % Doppler Steering Vector. 
     ta = ones(N,1);
  
    wmdb30 = zeros(M*N,M1,length(Ks));
    
for n=1:length(Ks)
    
     N1 = N - Ks(n) + 1;                         % Number of possible sub-appertures.
     g0 = exp(-1i*2*pi*fspt*(0:N1-1)).';
     g30 = chebwin(N1,30).*g0;           % 30-dB spatial taper.
     W30 = zeros(M*Ks(n),M1);
     GTSN = zeros(N,Ks(n));
     
     for p=0:Ks(n)-1
          Jp = [zeros(p,N1); eye(N1); zeros(N-N1-p,N1)];    % Beam Selector Matrix Jsm.
          GTSN(p+1:N1+p,p+1) = (Jp.'*Phi_jn*Jp)\g30;
     end
     
     utab30 = kron(tb.*bt,GTSN'*(ta.*at));         
     
     % W matrix construction.
     for p=0:M1-1
          Jp = [zeros(p,Kt); eye(Kt); zeros(M-Kt-p,Kt)];                  % Selector Matrix Jp.
          Rupdb30 = kron(Jp,GTSN)'*Ru*kron(Jp,GTSN);
          wp30 = Rupdb30\utab30;
          Wp30 = reshape(wp30,Ks(n),Kt);
          aux2 = Wp30*Jp.';
          W30(:,p+1)  = aux2(:);
     end   
   
     for m=1:M1
          wmdb30(:,m,n) = kron(eye(M),GTSN)*W30*F30(:,m);
     end

end

% LSINR Computations for Displaced Beams. 
SINR30_cube = zeros(M1,Lfd,length(Ks));
SINR30 = zeros(length(Ks),Lfd);
for n=1:length(Ks)
     for m=1:M1
          for n1=1:Lfd
               bt1 = exp(1i*2*pi*omegad(n1)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
               vt = kron(bt1,at);
               SINR30_cube(m,n1,n) = abs(wmdb30(:,m,n)'*vt)^2/real(wmdb30(:,m,n)'*Ru*wmdb30(:,m,n));
          end
     end
     SINR30(n,:) = max(abs(SINR30_cube(:,:,n)));
end    
    
LSINRpreTSN   = SINR30/SNRo;
    

%% 5. Beamspace post-Doppler: Displaced Filters with Two Step Nulling:

% Create Doppler Filterbank Matrix F for Displaced Filter Beamspace post-Doppler STAP:
dopplerfilterbank = linspace(0,300,M+1);
omegadopplerbank = dopplerfilterbank/fr;
F0 = zeros(M1,M);
for m=1:M
      F0(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));
end
   
F30 = diag(chebwin(M1,30))*F0;
wTSN = zeros(M*N,length(Ks));
  
for n=1:length(Ks)
      
     N1 = N - Ks(n) +1;                          % Number of possible sub-appertures.
        
% Beamformer G Matrix Construction for Displaced Filter Beamspace post-Doppler STAP
% with Two Step Nulling (TSN):
      taper_level = 30;                                            % Tapering level in db.
      g0 = exp(-1i*2*pi*fspt*(0:N1-1)).';
      gtap = chebwin(N1,taper_level).*g0;           % spatial taper application.

 % Calculation of Beamformer Matrix for Two-Step Nulling Approach:
      GTSN = zeros(N,Ks(n));
      for p=0:Ks(n)-1
           Jp = [zeros(p,N1); eye(N1); zeros(N-N1-p,N1)];    % Beam Selector Matrix Jsm.
           GTSN(p+1:N1+p,p+1) = (Jp.'*Phi_jn*Jp)\gtap;
      end

      SINRTSNdf_cube = zeros(M1,Lfd,length(Ks));

 % SINR Computations for Displaced Filter Beamspace post-Doppler STAP:
 % Solve a separate adaptive problem in each Doppler bin m:
      for m=1:M
           Fm30 = toeplitz([F30(:,m); zeros(Kt-1,1)],[F30(1,m) zeros(1,Kt-1)]);
           TmTSN = kron(Fm30,GTSN);
           Rum30 = TmTSN'*Ru*TmTSN;
           bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).'; 
           gt = kron(bdfb,at);
           utmTSM = TmTSN'*gt;
           wmTSN = Rum30\utmTSM;
           wTSN(:,n) = TmTSN*wmTSN;

           for n1=1:Lfd
                bt = exp(1i*2*pi*omegad(n1)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
                vt = kron(bt,at);
               SINRTSNdf_cube(m,n1,n) = abs(wTSN(:,n)'*vt)^2/real(wTSN(:,n)'*Ru*wTSN(:,n));
           end
       
      end
      
      SINR30(n,:) = max(abs(SINRTSNdf_cube(:,:,n)));
      
end

LSINRpostTSN = SINR30/SNRo;


%% Plot the SINR Losses:
     figure('NumberTitle', 'off','Name', ... 
                'Figure 67(a) SINR loss for example system, â=2.6, ó_õ = 0.2 m/sec, 10 deg vel. misalinment');
     plot(fd,10*log10(LSINRopt),'k','LineWidth',1.5); hold on;
     plot(fd,10*log10(LSINRsub2),'m','LineWidth',1.5);
     plot(fd,10*log10(LSINR30),'g','LineWidth',1.5);
     plot(fd,10*log10(LSINRpreTSN),'b','LineWidth',1.5);
     plot(fd,10*log10(LSINRpostTSN),'r','LineWidth',1.5);
     
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 1]);
     xlim([-5 305]);
     legend('Optimum Fully Adaptive','Elem-Pre (54)','Elem-Post (54)',...
                  'Beam-Pre TSN (15)','Beam-Post TSN (15)','Location','SouthEast')
     grid on;