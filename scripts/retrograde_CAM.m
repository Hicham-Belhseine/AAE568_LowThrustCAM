%% Parameters
% Timing
dtheta_t = 200 * pi/180;
n_revs = 2:.1:7;

dtheta_c = 2*pi*n_revs - dtheta_t/2;

% Vehicle
s1 = 8;  % Radius of vehicle 1 (m)
s2 = 8;  % Radius of vehicle 2 (m)

sigx = 80;    % x-position covariance
sigy = 1500;  % y-position covariance

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

T1 = 2*pi*sqrt(r1^3 / mu);
T2 = 2*pi*sqrt(r2^3 / mu);

n1 = 2*pi / T1;
n2 = 2*pi / T2;

xi_tf = 200 - 2*a0*(r1^3)/mu * ...
    (dtheta_t + sin(dtheta_c) - sin(dtheta_t+dtheta_c))

zeta_tf = a0*(r1^3)/mu*cos(kappa/2)* ...
    (3*dtheta_t*(dtheta_t/2 + dtheta_c)  -  8*sin(dtheta_t/2)*sin(dtheta_t/2 + dtheta_c))

plot(zeta_tf,xi_tf); grid on; box on;

bplane2x = @ (xi,zeta,Theta) xi*cos(Theta) + zeta*sin(Theta);
bplane2y = @ (xi,zeta,Theta) xi*sin(Theta) - zeta*cos(Theta);
x = bplane2x(xi_tf,zeta_tf,5*pi/180);
y = bplane2y(xi_tf,zeta_tf,5*pi/180);

semilogy(n_revs,p_collision(x,y,sigx,sigy,16))