%% Figure 51. SINR Loss performance for beamspace pre-Doppler STAP in a clutter-only
%% scenario.
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only.
% All rights reserved.

clc; clear; close all;

%% Radar System Ope3rational Parameters:
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
beta = 1;                        % beta parameter.
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
      ac(:,k) = exp(1i*2*pi*fspc(k)*(0:N-1));              % Clutter Spatial Steering Vector.
      bc(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1));      % Clutter Temporal Steering Vector.   
     Vc(:,k) = kron(bc(:,k),ac(:,k));                             % Clutter Space-Time Steering Vector.
end
     
Rc = real(Vc*Ksic*Vc');              
Rn = sigma2*eye(M*N);         
      
     %% Analytic Interference Covariance Matrix for Fully Adaptive STAP:
     Ru = Rc + Rn;    
     InvRu = inv(Ru);
     
  %% Target Space-Time Steering Vector:
     phit = 0; thetat = 0;                              % Target azimuth and elevation angles in degrees.
     fdt = 100;                                              % Target Doppler Frequency.
     fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
     omegat = fdt/fr;
     at = exp(-1i*2*pi*fspt*(0:N-1)).';           % Target Spatial Steering Vector.
     fd = 0:.5:300;   Lfd = length(fd);
     omegad = fd/fr;
     SNRo = M*N;

 %% LSINR Computation for Optimum Fully Adaptive Case:
   LSINRopt = zeros(1,Lfd);
     for n=1:Lfd
          bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Target Doppler Steering Vector   
          vt = kron(bt,at);
          w = InvRu*vt; %#ok<MINV>
          LSINRopt(n) = real(w'*vt)/SNRo;
     end
    
    %% Displaced Beam pre-Doppler Calculations:
    Kt =  2;                                        % sub-CPI length (fixed). 
    Ks = 2:6;                                     % various sub-apperture lengths or number of beams used. 
    M1 = M - Kt +1;                          % Number of possible sub-CPI's.

    %% Create Doppler Filterbank Matrix F:
     dopplerfilterbank = linspace(0,300,M1+1);
     omegadopplerbank = dopplerfilterbank/fr;
     F0 = zeros(M1,M1);
     for m1=1:M1
           F0(:,m1) = 1/sqrt(M1)*exp(-1i*2*pi*omegadopplerbank(m1)*(0:M1-1));
     end
   
     F30 = diag(chebwin(M1,30))*F0;
     
   %% Beamformer G Matrix Construction for Displaced Beam pre-Doppler.
     fsptd = 0;
     tb = [1;   1;];                                           % Temporal Binomial Taper.
     bt = [1;  -1;];                                           % Doppler Steering Vector. 
     ta = ones(N,1);
     at = ones(N,1);                                    % Steering Vector Pointing Broadside (phi = 0);
    wmdb0 = zeros(M*N,M1,length(Ks));
    wmdb30 = zeros(M*N,M1,length(Ks));
    
for n=1:length(Ks)
    
     N1 = N - Ks(n) + 1;                         % Number of possible sub-appertures.
     g0 = exp(-1i*2*pi*fsptd*(0:N1-1)).';
     g30 = chebwin(N1,30).*g0;           % 30-dB spatial taper.
     G0   = toeplitz([g0; zeros(Ks(n)-1,1);], [g0(1) zeros(1,Ks(n)-1)]);       % N1 x Ks Beamformer matrix G.  
     G30 = toeplitz([g30; zeros(Ks(n)-1,1);], [g30(1) zeros(1,Ks(n)-1)]);  % N1 x Ks Beamformer matrix G.  
     utab0 = kron(tb.*bt,G0'*(ta.*at));                                                            % Desired response. Eq. 236.
     utab30 = kron(tb.*bt,G30'*(ta.*at));                                                        % Desired response. Eq. 236.
     W0 = zeros(M*Ks(n),M1);
     W30 = zeros(M*Ks(n),M1);
          
     %% W matrix construction.
     for p=0:M1-1
          Jp = [zeros(p,Kt); eye(Kt); zeros(M-Kt-p,Kt)];                  % Selector Matrix Jp.
          Rupdb0   = kron(Jp,G0)'*Ru*kron(Jp,G0);
          Rupdb30 = kron(Jp,G30)'*Ru*kron(Jp,G30);
          wp0   = Rupdb0\utab0;
          wp30 = Rupdb30\utab30;
          Wp0   = reshape(wp0,Ks(n),Kt);
          Wp30 = reshape(wp30,Ks(n),Kt);
          aux1 = Wp0*Jp.';
          aux2 = Wp30*Jp.';
          W0(:,p+1)   = aux1(:);
          W30(:,p+1) = aux2(:);
     end   
   
     % Next calculate the Beamspace pre-Doppler composite vector for every Doppler bin m:
     for m=1:M1
          wmdb0(:,m,n)   = kron(eye(M),G0)*W0*F0(:,m);                % Eq. 246 for Untapered case.
          wmdb30(:,m,n) = kron(eye(M),G30)*W30*F30(:,m);         % Eq. 246 for Tapered case.
     end

end

 %% LSINR Computations for Displaced Beams. 
SINR0_cube   = zeros(M1,Lfd,length(Ks));
SINR30_cube = zeros(M1,Lfd,length(Ks));
SINR0 = zeros(length(Ks),Lfd);
SINR30 = zeros(length(Ks),Lfd);
for n=1:length(Ks)
     for m=1:M1
          for n1=1:Lfd
               bt1 = exp(-1i*2*pi*omegad(n1)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
               vt = kron(bt1,at);
               SINR0_cube(m,n1,n)   = abs(wmdb0(:,m,n)'*vt)^2/real(wmdb0(:,m,n)'*Ru*wmdb0(:,m,n));
               SINR30_cube(m,n1,n) = abs(wmdb30(:,m,n)'*vt)^2/real(wmdb30(:,m,n)'*Ru*wmdb30(:,m,n));
          end
     end
        SINR0(n,:)   = max(abs(SINR0_cube(:,:,n)));
        SINR30(n,:) = max(abs(SINR30_cube(:,:,n)));
end    
    
    LSINR0     = SINR0/SNRo;
    LSINR30   = SINR30/SNRo;
    
 %% LSINR Calculations for Adjacent Beam pre-Doppler STAP:   
   selbeamdist = 6.35;                     % Selected Beam Angular Distance in degrees.
   beamangles1 = [-9:-1  0:8]*selbeamdist; thetabeam = theta;
  
 % beamangles1 = linspace(-4.5,4,N);  thetabeam = theta;      % Inc = .5.
 % beamangles1 = linspace(-9,8,N);     thetabeam = theta;      % Inc = 1.
 % beamangles1 = linspace(-18,16,N); thetabeam = theta;      % Inc = 2.
 % beamangles1 = linspace(-27,24,N); thetabeam = theta;      % Inc = 3.
 % beamangles1 = linspace(-36,32,N); thetabeam = theta;      % Inc = 4.
 % beamangles1 = linspace(-45,40,N); thetabeam = theta;      % Inc = 5.
 % beamangles1 = linspace(-54,48,N); thetabeam = theta;      % Inc = 6.
 % beamangles1 = linspace(-63,56,N); thetabeam = theta;      % Inc = 7.
 % beamangles1 = linspace(-72,64,N); thetabeam = theta;      % Inc = 8.
 % beamangles1 = linspace(-81,72,N); thetabeam = theta;      % Inc = 9.
 % beamangles1 = linspace(-90,80,N); thetabeam = theta;      % Inc = 10.
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
 
 %% Plot Beams:
%  angle_edge = 90;
%  phi2 = -angle_edge:.2:angle_edge;  Lphi2 = length(phi2);
%  Resp_odd = zeros(N,Lphi2);
%  Resp_even = zeros(N,Lphi2);
%  
%  for n=1:N
%       for k= 1:Lphi2
%            fspt = d/lambda*sin(phi2(k)*pi/180);
%            at = exp(-1i*2*pi*fspt*(0:N-1)).';
%            Resp_odd(n,k)  = Godd(:,n).'*at;
%            Resp_even(n,k) = Geven(:,n).'*at;
%       end
%  end
%  [Azimuth Element] = meshgrid(phi2,1:N);
%  figure()
%  subplot(1,2,1);
%  pcolor(Azimuth, Element, abs(Resp_odd));
%  plot(phi2, abs(Resp_odd));
%  xlim([-angle_edge angle_edge]);
%  grid on;
%  shading flat;
%  hold on;
%  line([0 0], [0 5], [4.5 4.5],'Color','k', 'LineWidth',2,'EraseMode','normal');
%  title('Response of Odd Beams')
%  subplot(1,2,2)
%  pcolor(Azimuth, Element, abs(Resp_even));
%  plot(phi2, abs(Resp_even));
%  xlim([-angle_edge angle_edge]);
%  grid on;
%  shading flat;
%  hold on;
%  line([0 0], [0 5], [4.5 4.5],'Color','k', 'LineWidth',2,'EraseMode','normal');
%  title('Response of Even Beams');
 
 %%    
     
 wmab0   = zeros(M*N,M1,length(Ks));
 wmab30 = zeros(M*N,M1,length(Ks));
    
for n=1:length(Ks)
    
     N1 = N - Ks(n) + 1;                                                                       % Number of possible sub-appertures.
     P = floor(Ks(n)/2);
     if mod(Ks(n),2)  % if Ks is odd
             Js = [zeros(N/2-P,Ks(n)); eye(Ks(n)); zeros(N/2-P-1,Ks(n))];    % Beam Selector Matrix Js.
             Gab0 = Godd*Js;                                                                          % N x Ks Beamformer matrix.
             Gab30 = diag(chebwin(N,30))*Gab0;
     else
             Js = [zeros(N/2-P,Ks(n)); eye(Ks(n)); zeros(N/2-P,Ks(n))];        % Beam Selector Matrix Js.    
             Gab0 = Geven*Js;                                
             Gab30 = diag(chebwin(N,30))*Gab0;% N x Ks Beamformer matrix.
     end
  
%% Selected Adjacent Beams Response.
%    Respab0 = zeros(Ks(n),Lphi2);  
%    Respab30 = zeros(Ks(n),Lphi2);  
%    for n1=1:Ks(n)
%         for k= 1:Lphi2
%               fspt1 = d/lambda*sin(phi2(k)*pi/180);
%               at1 = exp(-1i*2*pi*fspt1*(0:N-1)).';
%               Respab0(n1,k)  = Gab0(:,n1).'*at1;
%               Respab30(n1,k)  = Gab30(:,n1).'*at1;
%         end
%    %     polardb(phi2*pi/180, 10*log10(abs(Respab(n1,:))), -30);
%    %    hold on;
%    % Normalization:
%    Respab0(n1,:) = Respab0(n1,:)/max(abs(Respab0(n1,:))); 
%    Respab30(n1,:) = Respab30(n1,:)/max(abs(Respab30(n1,:))); 
%    end
%          
%  % figure()
%  plot(phi2, abs(Respab0));
%  hold on;
%  plot(phi2, abs(Respab30),'r');
%  xlim([-angle_edge angle_edge]);
%  grid on;     
          
     %% W matrix construction.
     ut2ab0   = kron(tb.*bt,Gab0'*(ta.*at));                           % Desired response. Eq. 236.
     ut2ab30 = kron(tb.*bt,Gab30'*(ta.*at));                        % Desired response. Eq. 236.
     Wab0 = zeros(M*Ks(n),M1);
     Wab30 = zeros(M*Ks(n),M1);
     
     for p=0:M1-1
          Jp = [zeros(p,Kt); eye(Kt); zeros(M-Kt-p,Kt)];         % Selector Matrix Jp.
          Rup0   = kron(Jp,Gab0)'*Ru*kron(Jp,Gab0);
          Rup30 = kron(Jp,Gab30)'*Ru*kron(Jp,Gab30);
          wpab0 = Rup0\ut2ab0;
          wpab30 = Rup30\ut2ab30;
          Wpab0 = reshape(wpab0,Ks(n),Kt);
          Wpab30 = reshape(wpab30,Ks(n),Kt);
          aux1 = Wpab0*Jp.';
          aux2 = Wpab30*Jp.';
          Wab0(:,p+1)    = aux1(:);
          Wab30(:,p+1)  = aux2(:);
     end 

     for m=1:M1
          wmab0(:,m,n)   = kron(eye(M),Gab0)*Wab0*F0(:,m);
          wmab30(:,m,n) = kron(eye(M),Gab30)*Wab30*F30(:,m);
     end

end

 %% LSINR Computations for Displaced Beams. 
SINR0ab_mat   = zeros(M1,Lfd,length(Ks));
SINR30ab_mat   = zeros(M1,Lfd,length(Ks));
SINR0ab = zeros(length(Ks),Lfd);
SINR30ab = zeros(length(Ks),Lfd);
for n=1:length(Ks)
     for m=1:M1
           for n1=1:Lfd
                 bt = exp(-1i*2*pi*omegad(n1)*(0:M-1)).'; % Dummy Target Doppler Steering Vector   
                 vt = kron(bt,at);
                 SINR0ab_mat(m,n1,n)   = abs(wmab0(:,m,n)'*vt)^2/real(wmab0(:,m,n)'*Ru*wmab0(:,m,n));
                 SINR30ab_mat(m,n1,n) = abs(wmab30(:,m,n)'*vt)^2/real(wmab30(:,m,n)'*Ru*wmab30(:,m,n));
           end
     end
     SINR0ab(n,:)   = max(abs(SINR0ab_mat(:,:,n)));
     SINR30ab(n,:) = max(abs(SINR30ab_mat(:,:,n)));
end    
    
LSINR0ab     = SINR0ab/SNRo;
LSINR30ab   = SINR30ab/SNRo;

         %% Plot the SINR Losses:
     figure('NumberTitle', 'off','Name', ... 
         'Figure 51. SINR Loss Performance for Beamspace pre-Doppler STAP in a Clutter-only Scenario');
     
    %  Displaced Beams Results:
       subplot(2,2,1);
%    plot(fd,10*log10(LSINRopt),'LineWidth',1.5)      
%    hold on;
     plot(fd,10*log10(LSINR0),'LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 0]);
     xlim([-5 305]);
     legend('2 Beams', '3 Beams', '4 Beams', '5 Beams', '6 Beams','Location','South');
     title('Displaced Beams, Untapered');
     
     subplot(2,2,3);
  %   plot(fd,10*log10(LSINRopt),'LineWidth',1.5)      
  %   hold on;
     plot(fd,10*log10(LSINR30),'LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 0]);
     xlim([-5 305]);
     legend('2 Beams', '3 Beams', '4 Beams', '5 Beams', '6 Beams','Location','South');
     title('Displaced Beams, 30 dB Chebyshev Taper');
     
     
     % Adjacent Beams Results:
     subplot(2,2,2);
%     plot(fd,10*log10(LSINRopt),'LineWidth',1.5)      
 %    hold on;
     plot(fd,10*log10(LSINR0ab),'LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 0]);
     xlim([-5 305]);
     legend('2 Beams', '3 Beams', '4 Beams', '5 Beams', '6 Beams','Location','South');
     title('Adjacent Beams, Untapered');
     
     subplot(2,2,4);
  %   plot(fd,10*log10(LSINRopt),'LineWidth',1.5)      
  %   hold on;
     plot(fd,10*log10(LSINR30ab),'LineWidth',1.5)
     grid on;
     ylabel('SINR Loss (dB)');
     xlabel('Target Doppler Frequency (Hz)');
     ylim([-30 0]);
     xlim([-5 305]);
     legend('2 Beams', '3 Beams', '4 Beams', '5 Beams', '6 Beams','Location','South');
      title('Adjacent Beams, 30 dB Chebyshev Taper');