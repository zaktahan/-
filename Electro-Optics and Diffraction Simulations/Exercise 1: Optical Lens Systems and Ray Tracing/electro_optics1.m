%% Optical and Lens Simulation - Ilai Zaidel
% This MATLAB script calculates transmission, focal lengths, and ray tracing
% for different lenses and wavelengths using BK7 and F2 glass models.
% Author: Ilai Zaidel

%% --------- Setup ---------
close all; 
clc; 
clear;

%% --------- Question 1.d: S-Polarization Transmission ---------
n0 = 1.774; 
c = 3e8; 
lambda = linspace(1173e-9, 1177e-9, 1000); 
d = 0.00044879895; % Thickness
theta_i = atan(n0); 

a_values = [1e26, 1e27, 1e28]; 

figure;
hold on;

for i = 1:length(a_values)
    a = a_values(i);
    n = n0 - a*(lambda.^2)/c^2;           % Refractive index
    theta_t = asin(sin(theta_i)./n);      % Snell's law
    gamma_s = -sin(theta_i - theta_t)./sin(theta_i + theta_t); % S-polarization
    R = gamma_s.^2;
    delta = d*4*pi.*n.*cos(theta_t)./lambda;
    trans_ratio = (1-R).^2 ./ ((1-R).^2 + 4*R.*(sin(delta./2).^2));
    
    plot(lambda, trans_ratio);
end

title('S-Polarization Transmission vs. $\lambda$','Interpreter','latex');
xlabel('$\lambda \ [\mu {\rm m}]$','Interpreter','latex');
ylabel('Transmission','Interpreter','latex');
legend('$a=10^{26}$', '$a=10^{27}$', '$a=10^{28}$','Interpreter','latex');
hold off;

%% --------- Question 2.B: BiConvex Lens Focal Length ---------
R1 = 24e-3; 
R2 = -36e-3; 
lambda = 0.4:0.01:1; 

n = n_BK7(lambda); 
P = (n-1)./R1 + (1-n)./R2; 
f = 1./P; 

figure;
plot(lambda, f);
title('Focal Length vs. $\lambda$','Interpreter','latex');
xlabel('$\lambda \ [\mu {\rm m}]$','Interpreter','latex');
ylabel('$f \ [{\rm m}]$','Interpreter','latex');

%% --------- Question 2.C: Ray Tracing ---------
lambda1 = 0.624; 
lambda2 = 0.512; 
R = 24e-3; 
L1 = 70e-3; 

% Lambda1 calculation
n_lambda1 = n_BK7(lambda1);
L2_1 = -(1*L1 + 0) / (2*(1-n_lambda1)/R*L1 + 1);

figure;
hold on;
title('Ray Tracing for $-\pi/5 \le \theta \le \pi/5$, $x=1$','Interpreter','latex');
xlabel('$Z \ [{\rm m}]$','Interpreter','latex');
ylabel('Ray Height $[{\rm m}]$','Interpreter','latex');

for theta = -pi/5 : pi/25 : pi/5
    [ray, z] = bk7RayTrace(1, theta, lambda1, R, L1, L2_1);
    plot(z, ray);
end
hold off;

% Lambda1 and Lambda2 for x range
L2_2 = -(1*L1 + 0) / (2*(1-n_BK7(lambda2))/R*L1 + 1);

figure;
hold on;
title('Ray Tracing for $-1\le x\le1$, $\theta=0$','Interpreter','latex');
xlabel('$Z \ [{\rm m}]$','Interpreter','latex');
ylabel('Ray Height $[{\rm m}]$','Interpreter','latex');

theta = 0;
for x = -1:0.1:1
    [ray1, z1] = bk7RayTrace(x, theta, lambda1, R, L1, L2_1);
    plot(z1, ray1, 'r');
    [ray2, z2] = bk7RayTrace(x, theta, lambda2, R, L1, L2_2);
    plot(z2, ray2, 'b');
end
legend('$\lambda_1=0.624\ [\mu {\rm m}]$', '$\lambda_2=0.512\ [\mu {\rm m}]$', 'Interpreter','latex');
hold off;

%% --------- Question 2.D: Achromatic Lens Focal Length ---------
d1 = 0.43e-3; 
d2 = 2.3e-3;

syms R3;

% Lambda1
nF2_l1 = n_F2(lambda1);
nBK7_l1 = n_BK7(lambda1);

M1 = [1,0; (1-nBK7_l1)/(nBK7_l1*R), 1/nBK7_l1];
M2 = [1,d1;0,1];
M3 = [1,0; (nBK7_l1-nF2_l1)/(-nF2_l1*R), nBK7_l1/nF2_l1];
M4 = [1,d2;0,1];
M5 = [1,0; (nF2_l1-1)/(-R3), nF2_l1];

Mtot1 = M5*M4*M3*M2*M1;
C1 = Mtot1(2,1);

% Lambda2
nF2_l2 = n_F2(lambda2);
nBK7_l2 = n_BK7(lambda2);

M1 = [1,0; (1-nBK7_l2)/(nBK7_l2*R), 1/nBK7_l2];
M2 = [1,d1;0,1];
M3 = [1,0; (nBK7_l2-nF2_l2)/(-nF2_l2*R), nBK7_l2/nF2_l2];
M4 = [1,d2;0,1];
M5 = [1,0; (nF2_l2-1)/(-R3), nF2_l2];

Mtot2 = M5*M4*M3*M2*M1;
C2 = Mtot2(2,1);

% Solve R3 and convert to numeric
R3_val = double(vpa(solve(C1-C2==0, R3),7)); % numeric array
disp(['Calculated R3: ', num2str(R3_val(1))]); % first solution

%% --------- Functions ---------
function nBK7 = n_BK7(lambda)
% Refractive index for BK7 glass
nBK7Squared = 1 + 1.03961212.*(lambda.^2)./(lambda.^2 - 0.00600069867) ...
                 + 0.231792344.*(lambda.^2)./(lambda.^2 - 0.0200179144) ...
                 + 1.01046945.*(lambda.^2)./(lambda.^2 - 103.560653);
nBK7 = sqrt(nBK7Squared);
end

function nF2 = n_F2(lambda)
% Refractive index for F2 glass
nF2Squared = 1 + 1.34533359.*(lambda.^2)./(lambda.^2 - 0.00997743871) ...
                 + 0.209073176.*(lambda.^2)./(lambda.^2 - 0.0470450767) ...
                 + 0.937357162.*(lambda.^2)./(lambda.^2 - 111.886764);
nF2 = sqrt(nF2Squared);
end

function [ray, Z] = bk7RayTrace(x, theta, lambda, R, L1, L2)
% Ray tracing through a BK7 bi-convex lens
n = n_BK7(lambda);
c = 2*(1-n)/R;
M_lens = [1, 0; c, 1];
M_free = [1, L1; 0, 1];

z1 = 0:1e-6:L1;
z2 = L1:1e-6:L1+L2;
Z = [z1, z2];

ray_before = theta*z1 + x;
ray_after = zeros(size(z2));
ray_after(1) = ray_before(end);

for i = 2:length(z2)
    M_temp = [1, z2(i)-L1; 0,1];
    Mtot = M_temp*M_lens*M_free;
    ray_after(i) = Mtot(1,1)*x + Mtot(1,2)*theta;
end

ray = [ray_before, ray_after];
end
