%%% Calculation of the probe diffraction in a plasma - Waterfall


hollow = 0; % Choose 1 for a simulation with hollow plasma
            % 0 for a cylindrical plasma

waterfall = 1; % Choose 1 to run a waterfall (superpostition of different propagations)
               % 0 just to make a test for one z
               
plotting = 0; % Choose 1 to plot the intensity for each z during waterfall


    %Properties

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
alpha = pi/2;

%Plasma size
    % Channel radius (plasma radius if cylindrical plasma)
r0=200E-6; % !! increase a or reduce calib if you decrease r0 !!
    % For hollow plasma (ring width)
r1=5E-6; % !! increase a or reduce calib if you decrease r1 !!

%Accuracy 
a = 3000; % step
calib = 3E-1; % to calibrate delta_xi and delta_eta = L/a = z*calib/a


% Test for one z
if (waterfall==0)
    z = 0.3;
    %L=1E-2;
    L = z * calib;
    L_x=2*pi*lambda*a*z/L; %axes of the camera plane
    if hollow==1
        A=hollowPlasma_aperture(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a);
    else
        A=cylindrical_aperture(nLi,nPlasma,lambda,r0,z,L,alpha,a);
    end
%     A=A';
%     plot(atan(imag(A(:))./real(A(:))));
    U=fft2(A);
    U=fftshift(U);
    U=abs(U).^2;
    imagesc([-L_x L_x],[-L_x L_x],U); axis xy;
    colormap hot;
    colorbar;
end %if waterfall 0


% FFT Calculation for different z , plotting and waterfall

if (waterfall==1)
    
    %Propagation from -zmax m to +zmax m, with nbz iterations, with a step of 2*zmax/nbz
    zmax = 2;
    nbz = 100;
    step = 2*zmax/nbz;
    
    W = ones(a,nbz+1); %waterfall matrix
    
    for k=0:nbz
        
        %Propagation
        z=-zmax+step*k
        
        if (z<-0.2 || z>0.2)
            
            %Keep the same size of window for different z : L' = z'*L/z = z'*calib
            L=z*calib;
            
            if hollow==1
                A=hollowPlasma_aperture(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a);
            else
                A=cylindrical_aperture(nLi,nPlasma,lambda,r0,z,L,alpha,a);
            end
            
            % abs(A(225,245))
            % sum((abs(A(:))-1).^2)
            % plot(atan(imag(A(:))./real(A(:))));
            
            U=fft2(A)*z/(1i*lambda);
            U=fftshift(U);
            U=abs(U).^2;
            
            %Proportional factor in fresnel diffraction
            
            % z=linspace(-5,5,a);
            % y=linspace(-5,5,a);
            % for k =1:a
            %     for j=1:a
            %         U(k,j)=U(k,j)*1/(i*lambda*r0)*exp(2*i*pi/lambda*(r0+((z(j)^2+y(k)^2)/(2*r0))));
            %     end
            % end
            
            
            %Plotting
            
            if (plotting == 1)
                figure;
                L_x=abs(z)*calib;
                imagesc([-L_x L_x],[-L_x L_x],U); axis xy;
                % caxis([2 1E8]);
                colormap hot;
                colorbar;
            end
            
            %Waterfall
            W(:,k+1)=U(:,500);
            
        end %if
    end %for k
    
    L_x=2*pi*lambda*a*calib; %axes
    figure;
    pcolor(linspace(-zmax,zmax,101),linspace(-L_x,L_x,a),W); axis xy;
    shading interp;
    colormap hot;
    
end %if waterfall 1