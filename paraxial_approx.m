%Probe propagation in the paraxial approximation near z=0
%Helmholtz equation become: laplacien_transverse(A)+2*pi*k*dA/dz=0
% k = 2*pi*n(milieu)/lambda

%function small_propagation = paraxial_approx(nLi,nPlasma,lambda,r0,r1,h,L,alpha,a)

   %Properties

%Density
Np = 1E17 ; %Plasma
NLi = 1E10; %Lithium
%Wavelength
lambda = 800E-9;
lambda_plasma = 3.34E4/sqrt(Np);
lambda_Li = 670E-9;
%Index of refraction
nPlasma = sqrt(1 - lambda^2/lambda_plasma^2);
re=0.2;
nLi = 1 + (NLi*re/(2*pi))*0.744/(1/lambda_Li^2 - 1/lambda^2);
%nPlasma=0.999967999488;
%nLi=1.0004;
%Angle between plasma (z axis) and probe direction
alpha = pi/10;
%Step (accuracy) 
a = 100;
%Plasma channel size
r0=50E-6;
r1=5E-6; %Plasma size
sigma = 1E-2/8;


%Choose 1 for a hollow plasma simulation
hollow=0;
plotting=1;


%Definition de la grille
L=1E-3;
L_xi=L/2;
L_eta=L/2;

xi=linspace(-L_xi,L_xi,a);
xi=(meshgrid(xi));
eta=linspace(-L_eta,L_eta,a);
eta=(meshgrid(eta))';

dxi = 2*L_xi/a;
deta = 2*L_eta/a;

%Faisceau gaussien at z=0
R=xi.^2+eta.^2;
A=ones(a);
A=A.*exp(-R/(2*sigma^2));
imagesc([-L_xi L_xi],[-L_xi L_xi],abs(A).^2);title ('plop');axis xy;colorbar;

for k=0:10
    
    h=k*0.001
    
    %Redefinition de la grille
    L=(1+h)*L
    L_xi=L/2;
    L_eta=L/2;
    
    dxi = 2*L_xi/a;
    deta = 2*L_eta/a;
    
    %Calcul du Laplacien transverse
    
    %d^2(A)/dx^2
    d2A_dxi2=sparse(a,a);
    for k=2:(a-1)
        d2A_dxi2(:,k)=(A(:,k+1)-2*A(:,k)+A(:,k-1))/(dxi^2);
    end
    
    %d^2(A)/dy^2
    d2A_deta2=sparse(a,a);
    for k=2:(a-1)
        d2A_deta2(:,k)=(A(k+1,:)-2*A(k,:)+A(k-1,:))/(deta^2);
    end
    
    %dA/dz
    dA_dz=sparse(a,a);
    d1=r0^2-eta.^2;
    dA_dz=1i*lambda/(4*pi*nPlasma)*(d2A_dxi2+d2A_deta2);
    
    if (hollow==1)
        d2=r1^2-eta.^2;
        dA_dz(d1<0)=1i*lambda/(4*pi*nLi)*(d2A_dxi2(d1<0)+d2A_deta2(d1<0));
        dA_dz(d2>0)=1i*lambda/(4*pi*nLi)*(d2A_dxi2(d2>0)+d2A_deta2(d2>0));
        
    else
        dA_dz(d1<0)=1i*lambda/(4*pi*nLi)*(d2A_dxi2(d1<0)+d2A_deta2(d1<0));
        
    end %if hollow
    
    %Calculation of A(h) : propagation of h in plasma
    %A(z+h)=A(z)+h*dA/dz+O(h)
    
    A=(A+h*dA_dz);
    
    %Plotting
    if (plotting==1)
        U=abs(A).^2;
        figure;
        imagesc([-L_xi L_xi],[-L_xi L_xi],U); axis xy;
        colormap jet;
        colorbar;
    end %if plotting
    
end %for k
    