close all;
clc;
clear;

%% --------- Question 1: Circular Aperture ---------
R = 0.05;         % Radius of circle [m]
L = 0.2;          % Observation window size [m]
N = 200;           % Number of samples
dX = L/N;          % Sample spacing
dY = L/N;

% x and y axes
x_n = -L/2:dX:(L/2)-dX;
y_n = (-L/2:dY:(L/2)-dY).';

% Generate circular aperture
circle = zeros(N,N);
for i = 1:N
    for j = 1:N
        r_ij = sqrt(x_n(i)^2 + y_n(j)^2);
        circle(i,j) = circ(r_ij, R);
    end
end

% Plot circular aperture
figure;
imagesc(x_n, y_n, circle);
axis image;
axis xy;
title('Circular Aperture circ(r,R)','interpreter','latex');
xlabel('x [m]','interpreter','latex');
ylabel('y [m]','interpreter','latex');

%% --------- Question 1.D: FFT of Circular Aperture ---------
fx = -1/(2*dX):1/L:1/(2*dX)-1/L;
fy = fx;

circle_fft = fftshift(abs(fft2(circle)));

% Plot FFT using imagesc
figure;
imagesc(fx, fy, circle_fft);
axis image;
axis xy;
colorbar;
title('FFT of circ(r,R) using imagesc','interpreter','latex');
xlabel('$f_x$ [cycles/m]','interpreter','latex');
ylabel('$f_y$ [cycles/m]','interpreter','latex');

% Plot FFT using surf
figure;
surf(fx, fy, circle_fft);
axis xy;
camlight left;
lighting phong;
shading interp;
title(' FFT of circ(r,R) using surf','interpreter','latex');
xlabel('$f_x$ [cycles/m]','interpreter','latex');
ylabel('$f_y$ [cycles/m]','interpreter','latex');

%% --------- Question 2: Fraunhofer Diffraction ---------
lambda = 1007e-9;         % Wavelength [m]
z0 = 5*R^2/lambda;        % Observation distance

x = fx*lambda*z0;
[X, Y] = meshgrid(x, x);
r = sqrt(X.^2 + Y.^2);
G = R^2 .* jinc(R .* r / (lambda*z0));

% Field intensity
field_intensity = abs(G/(lambda*z0)).^2;

figure;
imagesc(x, x, field_intensity);
axis image; axis xy;
colorbar;
title(' Fraunhofer diffraction, field intensity at z_0','interpreter','latex');
xlabel('x [m]','interpreter','latex');
ylabel('y [m]','interpreter','latex');

%% --------- Question 2.D: Field intensity along y=0 ---------
[FWHM, x_interp, field_cut_interp] = interp_FWHM(field_intensity, x);

figure;
plot(x_interp, field_cut_interp);
hold on;
width = max(field_cut_interp)/2;
fwhm_line = width * ones(1,length(x_interp));
plot(x_interp, fwhm_line,'r');
hold off;
title(' Field Cut and FWHM at y=0','interpreter','latex');
xlabel('x [m]','interpreter','latex');
ylabel('Intensity','interpreter','latex');
legend('Field Cut','FWHM');

disp(['FWHM: ', num2str(FWHM)]);

%% --------- Question 3: Fresnel Diffraction ---------
[u2, x_prop] = propFresnel(circle,L,lambda,z0/50);
fresnel_intensity = abs(u2).^2;

figure;
imagesc(x_prop, x_prop, fresnel_intensity);
colorbar;
axis image; axis xy;
title(' Fresnel Diffraction at z=z_0/50','interpreter','latex');
xlabel('x [m]','interpreter','latex');
ylabel('y [m]','interpreter','latex');

[FWHM, x_interp, fresnel_cut_interp, index1, index2] = interp_FWHM(fresnel_intensity,x_prop);

figure;
plot(x_interp, fresnel_cut_interp);
hold on;
width = max(fresnel_cut_interp)/2;
plot(x_interp(index1:index2), width*ones(1,length(index1:index2)),'r');
hold off;
title(' Fresnel diffraction, field cut y=0','interpreter','latex');
xlabel('x [m]','interpreter','latex');
ylabel('Intensity','interpreter','latex');
legend('Field Cut','FWHM');
disp(['FWHM Fresnel z0/50: ', num2str(FWHM)]);

%% --------- Multiple distances ---------
z = [0.01*z0,0.1*z0,0.5*z0,z0,2*z0,10*z0];
multy = [0.01,0.1,0.5,1,2,10];
fwhm_vec = zeros(size(z));

figure;
for i = 1:length(z)
    subplot(2,6,i);
    [u2, x_prop] = propFresnel(circle,L,lambda,z(i));
    fresnel_intensity = abs(u2).^2;
    imagesc(x_prop,x_prop,fresnel_intensity);
    axis image; axis xy;
    title(['$z=$', num2str(multy(i)),'$\cdot z_0$'],'interpreter','latex');
    xlabel('X [m]'); ylabel('Y [m]');
    
    subplot(2,6,i+6);
    [FWHM, x_interp, fresnel_cut_interp] = interp_FWHM(fresnel_intensity,x_prop);
    plot(x_interp,fresnel_cut_interp);
    title(['$z=$', num2str(multy(i)),'$\cdot z_0$ for y=0'],'interpreter','latex');
    xlabel('X [m]'); ylabel('Intensity');
    fwhm_vec(i) = FWHM;
end

figure;
plot(z, fwhm_vec);
title(' FWHM vs. distance z','interpreter','latex');
xlabel('z [m]','interpreter','latex');
ylabel('FWHM [m]','interpreter','latex');
