%Probe propagation in the paraxial approximation near z=0
%Helmholtz equation become: laplacien_transverse(A)+2*pi*k*dA/dz=0
% k = 2*pi*n(milieu)/lambda

%function small_propagation = paraxial_approx(nLi,nPlasma,lambda,r0,r1,h,L,alpha,a)

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
alpha = pi/2;

%Plasma size
    % Channel radius (plasma radius if cylindrical plasma)
r0=50E-6; % !! increase a or reduce calib if you decrease r0 !!
    % For hollow plasma (ring width)
r1=25E-6; % !! increase a or reduce calib if you decrease r1 !!

%Accuracy 
a = 100; % step
calib = 5E-2; % to calibrate delta_xi and delta_eta = L/a = z*calib/a (4E-2/0.4)


%Choose 1 for a hollow plasma simulation
hollow=0;
plotting=0;
waterfall = 1;


%Definition de la grille
L=8E-3;
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
sigma = 1E-2/4;
A=ones(a);
A=A.*exp(-R/(2*sigma^2));
imagesc([-L_xi L_xi],[-L_xi L_xi],abs(A).^2);title ('plop');axis xy;colorbar;

W = ones(a,11);

for k=0:10
    
    h=k*0.01
    
    %Redefinition de la grille
    L= (1+h)*L;
    L_xi=L/2;
    L_eta=L/2;
    xi=linspace(-L_xi,L_xi,a);
    xi=(meshgrid(xi));
    eta=linspace(-L_eta,L_eta,a);
    eta=(meshgrid(eta))';
    
    dxi = 2*L_xi/a;
    deta = 2*L_eta/a;
    
    %Calcul du Laplacien transverse
    
    %d^2(A)/dx^2
    d2A_dxi2=sparse(a,a);
    for j=2:(a-1)
        d2A_dxi2(:,j)=(A(:,j+1)-2*A(:,j)+A(:,j-1))/(dxi^2);
    end
    
    %d^2(A)/dy^2
    d2A_deta2=sparse(a,a);
    for j=2:(a-1)
        d2A_deta2(:,j)=(A(j+1,:)-2*A(j,:)+A(j-1,:))/(deta^2);
    end
    
    status = h
    dA_dz=sparse(a,a);
    %d1=r0^2-eta.^2;
    %dA_dz=1i*lambda/(4*pi*nPlasma)*(d2A_dxi2+d2A_deta2);
    
    if (hollow==1)
        d2=r1^2-eta.^2;
        dA_dz(d1<0)=1i*lambda/(4*pi*nLi)*(d2A_dxi2(d1<0)+d2A_deta2(d1<0));
        dA_dz(d2>0)=1i*lambda/(4*pi*nLi)*(d2A_dxi2(d2>0)+d2A_deta2(d2>0));
        
    else
        
        ys=eta.^2;
        for j=1:a
            for i=1:a
                ys=eta.^2;
                d0=r0.^2-ys(i,j);
                d1=sqrt(d0)/(sin(alpha));
                dA_dz=1i*lambda/(4*pi*nLi)*(d2A_dxi2+d2A_deta2);
                if ( (d0>=0) )%&& (xi(i,j) <= d1-h) )
                    dA_dz=1i*lambda/(4*pi*nPlasma)*(d2A_dxi2+d2A_deta2);
                end %if
            end %for i
        end %for j
        
    end %if hollow
    
    %Calculation of A(h) : h propagation in plasma
    %A(z+h)=A(z)+h*dA/dz+O(h)
    
    A=(A+h*dA_dz);
    U=abs(A).^2;
    
    %Plotting
    if (plotting==1)
        
        figure;
        imagesc([-L_xi L_xi],[-L_xi L_xi],U); axis xy;
        colormap jet;
        colorbar;
    end %if plotting
    
    %Waterfall
    W(:,k+1)=U(:,a/2);
        
   
end %for k

if (waterfall==1)
figure;
L_x=2*pi*lambda*a*calib; %axes
pcolor(linspace(0,0.1,11),linspace(-L_x,L_x,a),W);axis xy;
shading interp;
colormap jet;
end

    