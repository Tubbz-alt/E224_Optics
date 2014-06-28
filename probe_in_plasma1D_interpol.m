%% Probe in plasma in 1D with interpolation (same size of window before propagation)
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
r0=100E-6;
    % For hollow plasma (ring width)
r1=25E-6; %

%Accuracy 
a = 100000; %number of points in the first window 
%%a=50.000, 24 points in plasma at z=2m
%%a=100.000, 50 points in plasma at z=2m
%%a=200.000, 100 points in plasma at z=2m
b = 1000; %number of points in the interpolation grid
calib = 5E-1; % to calibrate the size of the window before propagation

%% FFT calculation for different z locations

%Propagation from -zmax meters to +zmax meters, with nbz iterations, with a step of 2*zmax/nbz
zmax = 1;
nbz = 500;
step = 2*zmax/nbz;
z_min = 0.05;
L = 0.05*calib;
W_linear = ones(b,nbz+1); %waterfall matrix with linear interpolation
%W_spline = ones(1000,nbz+1); %waterfall matrix with spline interpolation

for k=0:nbz
    
    %Propagation
    z=-zmax+step*k
    
    if (z<-z_min || z>z_min)
        
        A=aperture_fastWaterfall(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a,hollow);
        
        U=fft(A);
        U=fftshift(U);
        U=sqrt(1/(z*lambda*1i*a))*U;
        U=abs(U).^2;
        
        % Interpolation
        L_y=2*pi*lambda*a*z/L; %new window after fft
        window=(linspace(L_y/2,-L_y/2,a))';
        grid=(linspace(0.08,-0.08,b))'; %the new grid within we do the interpolation
        U_linear = interp1(window,U,grid,'linear');
        %U_spline = interp1(window,U,grid,'spline');
        
        W_linear(:,k+1)=U_linear;
        %W_spline(:,k+1)=U_spline; 
        
    end %if
 end %for k
 
 
%% Plotting waterfall

Xaxis = linspace(-zmax,zmax,nbz+1);

figure;

pcolor(Xaxis,grid,W_linear);
shading interp;
colormap jet; 
colorbar; caxis([0.01E8 2.5E8]);
if (hollow==1) 
    title ({'Hollow Plasma' 'linear interpolation'},'FontSize', 12);
else
    title ({'Cylindrical Plasma' 'linear interpolation'},'FontSize', 12);
end

xlabel ('propagation in plasma z (m )', 'FontSize', 12);
ylabel ('Y ( m )', 'FontSize', 12);

leg1 = {'alpha = ' num2str(alpha,'%0.2f')};
leg2 = {'r0 = ' num2str(r0)};
leg3 = {'r1 = ' num2str(r1)};
leg4 = {'a = ' num2str(a)};
leg5 = {'nbz = ' num2str(nbz)};
leg6 = {'calib = ' num2str(calib)};
string = [leg1 leg2 leg3 leg4 leg5 leg6];
legend = text(zmax+zmax/3,-0.02,string);

set(legend,'backgroundcolor','w');

% figure;
% 
% pcolor(Xaxis,grid,W_spline);
% shading interp;
% colormap jet; 
% colorbar; %caxis([0.1E11 5E11]); %set(gca,'YLim',[-0.08 0.08]);
% if (hollow==1) 
%     title ({'Hollow Plasma' 'spline'},'FontSize', 12);
% else
%     title ({'Cylindrical Plasma' 'spline'},'FontSize', 12);
% end
% 
% 
% xlabel 'propagation in plasma z (m )';
% ylabel 'y ( m )';
% 
% leg1 = {'alpha = ' num2str(alpha,'%0.2f')};
% leg2 = {'r0 = ' num2str(r0)};
% leg3 = {'r1 = ' num2str(r1)};
% leg4 = {'a = ' num2str(a)};
% leg5 = {'nbz = ' num2str(nbz)};
% leg6 = {'calib = ' num2str(calib)};
% string = [leg1 leg2 leg3 leg4 leg5 leg6];
% legend = text(zmax+zmax/3,-0.05,string);
% 
% set(legend,'backgroundcolor','w');