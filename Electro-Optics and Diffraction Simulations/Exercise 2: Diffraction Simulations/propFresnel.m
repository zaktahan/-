function [u2, x_prop] = propFresnel(u1,L,lambda,z)
% Calculates Fresnel diffraction of 2D field u1 over distance z

N = size(u1,1);
dX = L/N;

x_tag = -L/2:dX:(L/2)-dX;
[X_tag, Y_tag] = meshgrid(x_tag);

k = 2*pi/lambda;

product = u1 .* exp(1i*k*(X_tag.^2 + Y_tag.^2)/(2*z));
fourier_product = fftshift(fft2(product));

fx = -1/(2*dX):1/L:1/(2*dX)-1/L;
x = fx * lambda * z;
[X,Y] = meshgrid(x);

x_prop = x;

u2 = exp(1i*k*z) .* exp(1i*k*(X.^2+Y.^2)/(2*z)) .* fourier_product / (1i*lambda*z);

end
