clc;clear;close all;

%% Parameters
% Vehicle
s1 = 8;  % Radius of vehicle 1 (m)
s2 = 1;  % Radius of vehicle 2 (m)

sigy_range = 50:5:300; % range for y positional covariance
sigx_range = sigy_range;

Theta = 5 * pi/180;  % covariance rotation wrt B-plane axis

m1 = 300;
f1 = 3e-3; 
a0 = f1 / m1;  % maximum acceleration magnitude

% Altitude/velocity
re = 6371e3;  % radius of Earth (m)

h1 = 550.0e3;  % Orbit height of vehicle 1 (m)
h2 = 550.2e3;  % Orbit height of vehicle 2 (m)

r1 = re + h1;  % Orbital radius of vehicle 1 (m)
r2 = re + h2;  % Orbital radius of vehicle 2 (m)

mu = 3.986e14; % Earth gravitational parameter (m^3/s^2)
v1 = sqrt(mu / r1);  % Orbital speed of vehicle 1 (m/s)

% Orbital params
i1 = 53 * pi/180;  % inclination of S1 (rad)
i2 = 53 * pi/180;  % inclination of S2 (rad)

RAAN1 = 180 * pi/180;  % RAAN of S1 (rad)
RAAN2 =   0 * pi/180;  % RAAN of S2 (rad)

% mutual inclination, angle between orbital planes
kappa = acos(sin(i1)*sin(i2)*cos(RAAN1-RAAN2) + cos(i1)*cos(i2));  

% Timing
T1 = 2*pi*sqrt(r1^3 / mu);
T2 = 2*pi*sqrt(r2^3 / mu);

n1 = 2*pi / T1;
n2 = 2*pi / T2;

dtheta_t = 20 * pi/180;
n_revs = 5;

% time = linspace(0,(2*pi*n_revs + dtheta_t/2)/n1,2e3);
time = linspace(0,dtheta_t/n1,5e2);
theta_c = 2*pi*n_revs + dtheta_t/2;
tc = time(end);

% Run TPBVP
solinit = bvpinit(time, [r2-r1, 0, 0.5, 0.5]);
options = bvpset('Stats','off','RelTol',1e-8);

p_tf = zeros(length(sigx_range),1);
p_t0 = zeros(length(sigx_range),1);
mean_sig = zeros(length(sigx_range),1);
for i = 1:length(sigx_range)
    sol = bvp4c( ...
        @ (t,y) CAM_ode(t,y,kappa,v1,n1,r1,a0,dtheta_t,theta_c), ...
        @ (ya,yb) CAM_bc(ya,yb,r2-r1,0,sigx_range(i),sigy_range(i),Theta), ...
        solinit, ...
        options  ...
    );
    
    xi = sol.y(1,:);
    zeta = sol.y(2,:);

    bplane2x = @ (xi,zeta,Theta) xi*cos(Theta) + zeta*sin(Theta);
    bplane2y = @ (xi,zeta,Theta) xi*sin(Theta) - zeta*cos(Theta);

    % bplane2x = @ (xi,zeta,Theta) zeta*cos(Theta) + xi*sin(Theta);
    % bplane2y = @ (xi,zeta,Theta) zeta*sin(Theta) - xi*cos(Theta);

    x = bplane2x(xi,zeta,Theta);
    y = bplane2y(xi,zeta,Theta);

    p_t0(i) = p_collision(x(1),y(1),sigx_range(i),sigy_range(i),s1+s2);
    p_tf(i) = p_collision(x(end),y(end),sigx_range(i),sigy_range(i),s1+s2);
    mean_sig(i) = sqrt(sigx_range(i)^2+sigy_range(i)^2);
end

%% Post-processing

figure(1)
semilogy(sigx_range,p_tf, 'linewidth', 2, 'color', 'red');
grid on;
title('Final Probability of Collision After Optimal CAM over \sigma Range');
xlabel('$\sigma$ [m]', 'Interpreter', 'latex');
ylabel('Final Probability of Collision after CAM');

improvement = zeros(length(sigx_range),1);
for i = 1:length(sigx_range)
    improvement(i) = 100*(p_t0(i)-p_tf(i))/p_t0(i);
end

figure(2)
plot(sigx_range,improvement, 'linewidth', 2, 'color', 'green');
grid on;
title('Probability of Collision Improvement After Optimal CAM over \sigma Range');
xlabel('$\sigma$ [m]', 'Interpreter', 'latex');
ylabel('Probability of Collision Imrpovement [%]');
%% Functions
function res = CAM_bc(ya,yb,xi0,zeta0,sig_x,sig_y,Theta)
%---------------------------------------------------------
% Description: Set of ODEs describing a low-thrust collision avoidance
% maneuver (CAM)
% 
% Inputs:
%     ya, yb - state initial and final conditions, respectively
% 
%     y = [
%         xi,
%         zeta,
%         lambda1,
%         lambda2,
%     ];
%---------------------------------------------------------
    sig_xi = sqrt( ...
        (sig_x^2 + sig_y^2 - (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
    );

    sig_zeta = sqrt( ...
        (sig_x^2 + sig_y^2 + (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
    );

    rho_xi_zeta = (sig_x^2 - sig_y^2) * sin(2*Theta) / sig_xi / sig_zeta;
    
    Q = [
        1/sig_xi^2    -rho_xi_zeta/sig_xi/sig_zeta;
        -rho_xi_zeta/sig_xi/sig_zeta  1/sig_zeta^2;
    ];

%     Q = eye(2);

    % final time condition from Pontryagin's maximum principle
    lambda_tf = Q*[yb(1);yb(2)];

    res = [
        ya(1) - xi0;
        ya(2) - zeta0;
        yb(3) - lambda_tf(1);
        yb(4) - lambda_tf(2);
    ];

end