function U = fresnel1D(lambda, transmission, propag, L_window, nbpoints)

%Grid definition
grid = (linspace(L_window/2,-L_window/2,nbpoints))';

%Fresnel term
U = transmission .* exp(1i*pi/(lambda*propag)*grid.^2);

%FFT calculation
U=fft(U);
U=fftshift(U);

%Outside-FFT term
new_length = 2*pi*lambda*nbpoints/L_window * propag;
new_grid = (linspace(new_length/2,-new_length/2,nbpoints))'; %Grid after FFT

U=sqrt(L_window/(propag*lambda*1i*nbpoints))*U;
U = U .* exp(2*1i*pi/lambda*(propag+(new_grid.^2)/(2*propag)));

end