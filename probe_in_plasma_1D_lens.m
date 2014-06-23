%%% Probe in plasma in 1D with the lens
%%% !! Just for fast waterfall : if ou want to plot intensity at each z, run
%%% probe_in_plasma !!

hollow = 0; % Choose 1 for a simulation with hollow plasma
            % 0 for a cylindrical plasma

            
%Properties

%Density
Np = 1E17 ; %Plasma
NLi = 3E16; %Lithium

%Wavelength
lambda = 800E-9;
lambda_plasma = 3.34E4/sqrt(Np);
lambda_Li = 670E-9;

%Refractive index
nPlasma = sqrt(1 - lambda^2/lambda_plasma^2); %nPlasma=0.999967999488;
re=1E-7;
nLi = 1 + (NLi*re/(2*pi))*0.744/(1/lambda_Li^2 - 1/lambda^2); %nLi=1.0004;

%Angle between plasma (z axis) and probe direction
alpha = pi/10;

%Plasma size
% Channel radius (plasma radius if cylindrical plasma)
r0=100E-6; % !! increase a or reduce calib if you decrease r0 !!
% For hollow plasma (ring width)
r1=25E-6; % !! increase a or reduce calib if you decrease r1 !!

%Accuracy
a = 30000; % step
%%a=50.000, 24 points in plasma at z=2m
%%a=100.000, 50 points in plasma at z=2m
%%a=200.000, 100 points in plasma at z=2m
calib = 1E-1; % to calibrate delta_xi and delta_eta = L/a = z*calib/a (4E-2/0.4)

%Lens properties
f=1; %focus
p=1E-1; %radius


%%FFT Calculation for different z and waterfall

%Propagation from -zmax m to +zmax m, with nbz iterations, with a step of 2*zmax/nbz
zmax = 1;
nbz = 100;
step = 2*zmax/nbz;

W = ones(a,nbz+1); %waterfall matrix
L_y=2*pi*lambda*a/calib; %axes after first FFT and before lens
%L_y=5E-2;
L_f= 2*f*calib ; %2*pi*lambda*a*f/L_y; %axes after the second FFT (after the lens)

for k=0:nbz
    
    %Propagation in plasma
    z=-zmax+step*k
    
    if (z<-0.1 || z>0.1)
        
        %Keep the same size of window for different z : L' = z'*L/z = z'*b
        L=abs(z)*calib;
        
        %%Fresnel with a 2-z propagation from plasma to lens
        
        % Grid
        L_eta=L/2;
        eta=(linspace(L_eta,-L_eta,a))';
        
        % Gaussian beam
        sigma = calib/4;
        transmission = exp(-(eta.^2)/(2*sigma^2));
        
        % Transmission function calculation
        
        if (hollow == 1) %hollow plasma
            
            d0=r0^2-eta.^2;
            d1=eta.^2-(r0-r1)^2;
            
            transmission = transmission.*exp(-2*1i*pi/lambda*2/sin(alpha)*nLi*r0);
            transmission(d0>=0) = transmission(d0>=0).*exp(-2*1i*pi/lambda*2/sin(alpha)*sqrt(d0(d0>=0))*(nPlasma-nLi));
            transmission(d1<0) = transmission(d1<0).*exp(-2*1i*pi/lambda*2/sin(alpha)*sqrt(-d1(d1<0))*(nLi-nPlasma));
            
            
        else %cylindrical plasma
            
            Lp = 2*r0/sin(alpha);
            d0=r0^2-eta.^2;
            
            transmission = transmission * (exp(-2*1i*pi/lambda*nLi*Lp));
            transmission(d0>=0) = transmission(d0>=0) .* exp(-2*1i*pi/lambda*2/abs(sin(alpha))*(nPlasma-nLi).*sqrt(d0(d0>=0)));
            
        end
        
        % Fresnel term
        transmission = transmission .* exp(1i*pi/(lambda*(2-z))*eta.^2);
        
        U=fft(transmission);
        U=fftshift(U);
        hold on; plot(abs(U));
        %Phase factor in fresnel diffraction outside fft
        U = U.* (sqrt((2-z)/(lambda*1i)) * exp(2*1i*pi/lambda*((2-z)+(eta.^2)/(2*(2-z)))));
     
        %%Lens phase shift
        y=(linspace(L_y,-L_y,a))';
        %U(abs(y)>p)=0; %lens clipping
        Ul = U.*exp(-1i*pi/(lambda*f)*y.^2); %lens term
        
        %%Fresnel with a 2f propagation from lens to camera
        Uf = Ul .* exp(1i*pi/(lambda*(2*f))*y.^2);
%                 Uf=fft(Uf)*sqrt(2*f/(lambda*1i));
%                 Uf=fftshift(Uf);
        Uf=abs(Uf).^2;
        
        % Waterfall
        W(:,k+1)=Uf;
        
        
     end %if
  end %for k

% figure;
% pcolor(linspace(-zmax,zmax,nbz+1),linspace(-L_f,L_f,a),W);axis xy;
% %caxis([0.1E10 9E10]);
% shading interp;
% colormap jet;
% %set(gca,'YLim',[-1E-3 1E-3]);
