%%% Transmission function * Fresnel exponential

function transmission = hollowPlasma_aperture(nLi,nPlasma,lambda,r0,r1,z,L,alpha,a)

%Definition of (xi,eta) plane

L_eta=L/2;
L_xi=L/2;

eta=linspace(-L_eta,L_eta,a);
eta=(meshgrid(eta));
xi=linspace(-L_xi,L_xi,a);
xi=(meshgrid(xi))';

%Faisceau gaussien
R=eta.^2+xi.^2;
sigma = 1E-2/4;
transmission = exp(-R/(2*sigma^2));
%imagesc(transmission)

for k=1:a
    
    y2=xi(k,1)^2;
    d0=r0^2-y2;
    d1=y2-(r0-r1)^2;
    
    if (d0<0)
        transmission(k,:)=transmission(k,:).*exp(-2*1i*pi/lambda*2/sin(alpha)*nLi*r0);
    elseif (d0*d1>=0)
        transmission(k,:)=transmission(k,:).*exp(-2*1i*pi/lambda*2/sin(alpha)*(nLi*(r0-sqrt(d0))+nPlasma*(sqrt(d0))));
    else
        transmission(k,:)=transmission(k,:).*exp(-2*1i*pi/lambda*2/sin(alpha)*(nLi*(r0-sqrt(d0)+sqrt(-d1))+nPlasma*(sqrt(d0)-sqrt(-d1))));
    end %if
    
end %for k

%Fresnel term

transmission = transmission .* exp(1i*pi/(lambda*z)*R);

end