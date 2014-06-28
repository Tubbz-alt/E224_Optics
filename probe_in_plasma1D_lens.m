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
alpha = pi/300;

%Plasma size
% Channel radius (plasma radius if cylindrical plasma)
r0=100E-6; % !! increase a or reduce calib if you decrease r0 !!
% For hollow plasma (ring width)
r1=25E-6; % !! increase a or reduce calib if you decrease r1 !!

%Accuracy 
a = 100000; %number of points in the first window
b = 1000; %number of points in the interpolation grid
calib = 5E-1; % to calibrate the size of the window before propagation

%Lens properties
f=1; %focus
p=1E-1; %radius


%% FFT Calculation for different z locations

%Propagation in plasma from -zmax meters to +zmax meters, with nbz iterations, with a step of 2*zmax/nbz
zmax = 1;
nbz = 500;
step = 2*zmax/nbz;
z_min = 0.05;

L_bp = 0.05*calib; %before plasma
eta = (linspace(L_bp,-L_bp,a))'; %first grid

W_linear = ones(b,nbz+1); %waterfall matrix with linear interpolation
%W_spline = ones(1000,nbz+1); %waterfall matrix with spline interpolation

for k=0:nbz
    
    %Propagation in plasma
    z=-zmax+step*k
    
    if (z<-z_min || z>z_min)
                 
        
        %% Fresnel in plasma
        
        %Gaussian beam
        sigma = 4E-2/16;
        transmission = exp(-(eta.^2)/(2*sigma^2));
        
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
        U = fresnel1D(lambda,transmission,z,L_bp,a);
       
        
        %% Fresnel from plasma to lens
        
        propag = 2*f; %-z;
        L_ap = 2*pi*lambda*a/L_bp * z;
        U = fresnel1D(lambda,U,propag,L_ap,a);

        %% Lens term
        
        L_lens = 2*pi*lambda*a/L_ap * propag;
        new_grid = (linspace(-L_lens/2,L_lens/2,a))';
        %U(abs(new_grid)>p) = 0; %lens clipping
        U = U.*exp(-1i*pi/(lambda*f)*new_grid.^2); %lens phase shift
        
        %% Fresnel from lens to camera
        
        U = fresnel1D(lambda,U,2*f,L_lens,a);
     
        U=abs(U).^2;

        % Interpolation
        L_cam = 2*pi*lambda*a*2*f/L_lens;
        window=(linspace(L_cam/2,-L_cam/2,a))';
        interp_grid=(linspace(0.04,-0.04,b))'; %the new grid within we do the interpolation
        U_linear = interp1(window,U,interp_grid,'linear');
        
        % Waterfall
        W_linear(:,k+1)=U_linear;
        
        
    end %if
end %for k


%% Plotting waterfall

figure;

pcolor(linspace(-zmax,zmax,nbz+1),interp_grid,W_linear); shading interp;
colorbar; %caxis([0.1E10 9E10]);
colormap jet;

