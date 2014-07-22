%% Probe in plasma in 1D with the lens
%%% !! Just for fast waterfall : if ou want to plot intensity at each z, run
%%% probe_in_plasma !!

hollow = 0; % Choose 1 for a simulation with hollow plasma
            % 0 for a cylindrical plasma


%% Properties

%Density
Np = 1E16 ; %Plasma
NLi = 1E16; %Lithium

%Wavelength
lambda = 800E-9;
lambda_plasma = 3.34E4/sqrt(Np);
lambda_Li = 670E-9;

%Refractive index
nPlasma = sqrt(1 - lambda^2/lambda_plasma^2); %nPlasma=0.999967999488;
re=2E-7;
nLi = 1 + (NLi*re/(2*pi))*0.744/(1/lambda_Li^2 - 1/lambda^2); %nLi=1.0004;

%Angle between plasma (z axis) and probe direction
alpha = pi/100;

%Plasma size
% Channel radius (plasma radius if cylindrical plasma)
r0=100E-6;
% For hollow plasma (ring width)
r1=25E-6;

%Accuracy
a = 100001; %number of points in the first window
b = 1000; %number of points in the interpolation grid
calib = 5E-1; % to calibrate the size of the window before propagation

%Lens properties
f=1; %focus
p=1.5E-2; %radius


%% FFT Calculation for different z locations

sigma = 1E-2/4; % RMS du faisceau gaussien

%Propagation in plasma from -zmax meters to +zmax meters, with nbz iterations, with a step of 2*zmax/nbz
zmax = 1;
nbz = 501;
step = 2*zmax/nbz;
L_bp = 0.08*calib; %before plasma
eta = (linspace(-L_bp/2,L_bp/2,a))'; %first grid

W_linear = ones(b,nbz+1); %waterfall matrix with linear interpolation

for k=0:nbz
    
    z=-zmax+step*k
    
    %% Diffraction in the plasma, propagation to the lens
    
    %Gaussian beam
    transmission = exp(-(eta.^2)/(2*(sqrt(2)*sigma)^2));
    %     figure(1);
    %     plot(eta,transmission.^2,'r');
    
    %Plasma transmission
    
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
    
    %Fresnel
    propag = 2*f-z;
    [U,new_grid] = fresnel1D(lambda,transmission,propag,eta);
    
    %% Lens term
    
    %lens clipping
    U(abs(new_grid)>p) = 0;
    %lens phase shift
    U = U.*exp(-1i*pi/(lambda*f)*new_grid.^2);
    
    %         hold on;
    %         %plot(new_grid,atan(Im(U)/r);
    %         plot(new_grid,abs(U).^2);
    %         xlim([-0.01 0.01]);
    %         pause(0.2);
    %         hold off;
    
    %% Fresnel propagation from lens to camera
    
    [U,grid_cam] = fresnel1D(lambda,U,2*f,new_grid);
    U=abs(U).^2;
    
    % Interpolation
    interp_grid=(linspace(-0.005,0.005,b))'; %the new grid within we do the interpolation
    U_linear = interp1(grid_cam,U,interp_grid,'linear');
    
    % Waterfall
    W_linear(:,k+1)=U_linear;
    
end %for k


%% Plotting waterfall

figure(2);

pcolor(linspace(-zmax,zmax,nbz+1),interp_grid,W_linear); shading interp;
axis xy;
colorbar; caxis([0.1 2]);
colormap jet;

