%%% Aperture in 1D for fast waterfall


function transmission = aperture_fastWaterfall(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a,hollow)
            
% y-profile of the probe at x=0 and before plasma

%Grid
L_y=L/2;
y=(linspace(L_y,-L_y,a))';

%Gaussian beam
sigma = 4E-2/16;
transmission = exp(-(y.^2)/(2*sigma^2));

%plot (transmission);

% Transmission function calculation

if (hollow == 1) %hollow plasma
    
    d0=r0^2-y.^2;
    d1=y.^2-(r0-r1)^2;
    
    transmission = transmission.*exp(-2*1i*pi/lambda*2/sin(alpha)*nLi*r0);
    transmission(d0>=0) = transmission(d0>=0).*exp(-2*1i*pi/lambda*2/sin(alpha)*sqrt(d0(d0>=0))*(nPlasma-nLi));
    transmission(d1<0) = transmission(d1<0).*exp(-2*1i*pi/lambda*2/sin(alpha)*sqrt(-d1(d1<0))*(nLi-nPlasma));
    
    
else %cylindrical plasma
    
    Lp = 2*r0/sin(alpha);
    d0=r0^2-y.^2;
    
    transmission = transmission * (exp(-2*1i*pi/lambda*nLi*Lp));
    transmission(d0>=0) = transmission(d0>=0) .* exp(-2*1i*pi/lambda*2/abs(sin(alpha))*(nPlasma-nLi).*sqrt(d0(d0>=0)));
    
end

%Fresnel term

transmission = transmission .* exp(1i*pi/(lambda*z)*y.^2);
 %plot(atan(imag(transmission)./real(transmission)));
% grid;
end