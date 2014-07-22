function [U,new_grid] = fresnel1D(lambda, transmission, propag, grid)

%Fresnel term
U = transmission .* exp(1i*pi/(lambda*propag)*grid.^2);

%FFT calculation
U=fft(U);
U=fftshift(U);

%New grid
new_grid = linspace(-lambda*propag/(grid(2)-grid(1))/2,lambda*propag/(grid(2)-grid(1))/2,size(U,1))';

%Outside-FFT term
U=(grid(2)-grid(1))/sqrt(propag*lambda*1i)*U;
U = U .* exp(2*1i*pi/lambda*(propag+(new_grid.^2)/(2*propag))) .*exp(-1i*2*pi/(lambda*propag)*grid(1)*new_grid);

end