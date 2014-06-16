%%% Transmission function * Fresnel exponential (neglected si Fraunhofer)

function transmission = cylindrical_aperture(nLi,nPlasma,lambda,r0,z,L,alpha,a)

%Definition of (xi,eta) plane

L_eta=L/2;
L_xi=L/2;

eta=linspace(-L_eta,L_eta,a);
eta=(meshgrid(eta));
xi=linspace(-L_xi,L_xi,a);
xi=(meshgrid(xi))';

y=xi.^2;
y=r0^2-y;
y(y<0)=0;

%Faisceau gaussien
% fen=5E-3;
% g1=meshgrid(linspace(-fen,fen,1000));
% g2=meshgrid(linspace(fen,-fen,1000))';
% G=g1.^2+g2.^2;
R=eta.^2+xi.^2;
sigma = 4E-2/32;

%Plasma (crossed by probe at y=0) length
Lp = 2*r0/sin(alpha);

transmission = exp(-R/(2*sigma^2)).*ones(a)*(exp(-2*1i*pi/lambda*nLi*Lp)); 
transmission = transmission .* exp(2*sqrt(y).*(-2*1i*pi/lambda*1/abs(sin(alpha))*(nPlasma-nLi)));
transmission = transmission .* exp(1i*pi/(lambda*z)*R);

%k=2*pi/(deltax)
%k=1/(lambda*z)

% for k=1:a   
%     
%     if ( ( (eta(k))>(-r0) ) && ( (eta(k))<r0 ) )
%         %transmission((a+1)-k,:) = transmission((a+1)-k,:)*exp(-2*i*pi/lambda*(1/abs(sin(alpha))*sqrt(r0^2-eta(k)^2)*(nPlasma-nLi)));
%         phase_shift((a+1)-k,:) = exp(-2*i*pi/lambda*(1/abs(sin(alpha))*sqrt(r0^2-eta(k)^2)*(nPlasma-nLi)));
%         %phase_shift((a+1)-k,:) = exp(-2*i*pi/lambda*nLi*L);
%         %phase_shift(k,:) = 1;
%         
%     %else
%         %phase_shift((a+1)-k,:) = exp(-2*i*pi/lambda*(1/abs(sin(alpha))*sqrt(r0^2-eta(k)^2)*(nPlasma-nLi)+L*nLi));
%         %phase_shift(k,:) = exp(-2*i*pi/lambda*(1/abs(sin(alpha))*sqrt(r0^2-eta(k)^2)*(nPlasma)));
%  
%     end
%     
%     for j=1:a
%         transmission((a+1)-k,j) = phase_shift((a+1)-k,j)*exp((i*pi/lambda)*(xi(j)^2+eta(k)^2)/r1);
%     
%     end
%     
% end

end
    