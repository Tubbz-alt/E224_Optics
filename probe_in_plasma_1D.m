%%% Probe in plasma in 1D
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
alpha = pi/300;

%Plasma size
    % Channel radius (plasma radius if cylindrical plasma)
r0=100E-6; % !! increase a or reduce calib if you decrease r0 !!
    % For hollow plasma (ring width)
r1=25E-6; % !! increase a or reduce calib if you decrease r1!!

%Accuracy 
a = 50000; % step 
%%a=50.000, 24 points in plasma at z=2m
%%a=100.000, 50 points in plasma at z=2m
%%a=200.000, 100 points in plasma at z=2m
calib = 1E-1; % to calibrate delta_xi and delta_eta = L/a = z*calib/a (4E-2/0.4)


%FFT Calculation for different z and waterfall

%Propagation from -zmax m to +zmax m, with nbz iterations, with a step of 2*zmax/nbz
zmax = 2;
nbz = 500;
step = 2*zmax/nbz;

W = ones(a,nbz+1); %waterfall matrix

for k=0:nbz
    
    %Propagation
    z=-zmax+step*k
    
    if (z<-0.1 || z>0.1)
        
        %Keep the same size of window for different z : L' = z'*L/z = z'*b
        L=abs(z)*calib;
        A=aperture_fastWaterfall(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a,hollow);
        
%         figure;
%         plot(abs(A));
%         
        % abs(A(225,245))
        % sum((abs(A(:))-1).^2)
        % plot(atan(imag(A(:))./real(A(:))));
        
        U=fft(A);
        U=fftshift(U);
        U=sqrt(z/(lambda*1i))*U;
        U=abs(U).^2;
       
        %Proportional factor in fresnel diffraction
        
%         for k =1:a
%             for j=1:a
%                 U(k,j)=U(k,j)*1/(i*lambda*r0)*exp(2*i*pi/lambda*(r0+((z(j)^2+y(k)^2)/(2*r0))));
%             end
%         end
        
        
        % Waterfall
        W(:,k+1)=U;
        
    end %if
 end %for k
 
figure;
L_x=2*pi*lambda*a/calib; %axes
pcolor(linspace(-zmax,zmax,nbz+1),linspace(-L_x,L_x,a),W);axis xy;
shading interp;
colormap jet; 
colorbar; caxis([0.1E11 5E11]);
if (hollow==1) 
    title 'hollow plasma'
else
    title 'cylindrical plasma'
end
xlabel 'propagation in plasma z (m )';
ylabel 'y ( m )';
set(gca,'YLim',[-0.08 0.08]);

% leg1 = {'nPlasma = ' num2str(nPlasma)};
% leg2 = {'nLi = ' num2str(nLi)};
leg3 = {'alpha = ' num2str(alpha,'%0.2f')};
leg4 = {'r0 = ' num2str(r0)};
leg5 = {'r1 = ' num2str(r1)};
leg6 = {'a = ' num2str(a)};
leg7 = {'nbz = ' num2str(nbz)};
leg8 = {'calib = ' num2str(calib)};
string = [leg3 leg4 leg5 leg6 leg7 leg8];
legend = text(2.6,-0.05,string);

set(legend,'backgroundcolor','w');











