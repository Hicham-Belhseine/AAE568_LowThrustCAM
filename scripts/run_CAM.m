%% Scenario Description
% Section VII Scenario in Hernando-Ayuso/Bombardelli paper

%% Parameters
% Vehicle
s1 = 8;  % Radius of vehicle 1 (m)
s2 = 8;  % Radius of vehicle 2 (m)

sig_x = 80;    % x-position covariance
sig_y = 1500;  % y-position covariance

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

dtheta_t = 200 * pi/180;
n_revs = 5;

% time = linspace(0,(2*pi*n_revs + dtheta_t/2)/n1,2e3);
time = linspace(0,dtheta_t/n1,5e2);
theta_c = 2*pi*n_revs + dtheta_t/2;
tc = time(end);

% Quadratic Cost Function
sig_xi = sqrt( ...
    (sig_x^2 + sig_y^2 + (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
);

sig_zeta = sqrt( ...
    (sig_x^2 + sig_y^2 - (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
);

rho_xi_zeta = (sig_x^2 - sig_y^2) * sin(2*Theta) / sig_xi / sig_zeta;

Q = [
    1/sig_xi^2    -rho_xi_zeta/sig_xi/sig_zeta;
    -rho_xi_zeta/sig_xi/sig_zeta  1/sig_zeta^2;
];

% Run TPBVP once
solinit = bvpinit(time, [r2-r1, 0, 1e-3, 1e-3]);
options = bvpset('Stats','off','RelTol',1e-8);

sol = bvp4c( ...
    @ (t,y) CAM_ode(t,y,kappa,v1,n1,r1,a0,dtheta_t,theta_c), ...
    @ (ya,yb) CAM_bc(ya,yb,r2-r1,0,Q), ...
    solinit, ...
    options  ...
);

t = sol.x*n1/2/pi;
xi = sol.y(1,:);
zeta = sol.y(2,:);

bplane2x = @ (xi,zeta,Theta) xi*cos(Theta) + zeta*sin(Theta);
bplane2y = @ (xi,zeta,Theta) xi*sin(Theta) - zeta*cos(Theta);

x = bplane2x(xi,zeta,Theta);
y = bplane2y(xi,zeta,Theta);

p_t0 = p_collision(x(1),y(1),sigx,sigy,s1+s2);
p_tf = p_collision(x(end),y(end),sigx,sigy,s1+s2);

revs = 0.5:.05:7;
p_optimal = zeros(1,length(revs));
p_retro = zeros(1,length(revs));
p_pro = zeros(1,length(revs));
for i = 1:length(revs)
    n_revs = revs(i);
    
    time = linspace(0,dtheta_t/n1,100);
    theta_c = 2*pi*n_revs + dtheta_t/2;

    % Costate guesses
    [lambda_retro_guess, lambda_pro_guess, state_retro_guess, state_pro_guess] = ...
        guess_CAM_costate(n_revs, dtheta_t, a0, r1, kappa, Q, r2-r1, 0);

    solinit_retro = bvpinit(time, [r2-r1, 0, lambda_retro_guess']);
    solinit_pro = bvpinit(time, [r2-r1, 0, lambda_pro_guess']);

    % Retrograde guess solution
    sol_retro = bvp4c( ...
        @ (t,y) CAM_ode(t,y,kappa,v1,n1,r1,a0,dtheta_t,theta_c), ...
        @ (ya,yb) CAM_bc(ya,yb,r2-r1,0,Q), ...
        solinit_retro, ...
        options  ...
    );

    t_retro = sol_retro.x*n1/2/pi;
    xi_retro = sol_retro.y(1,:);
    zeta_retro = sol_retro.y(2,:);

    x = bplane2x(xi_retro,zeta_retro,Theta);
    y = bplane2y(xi_retro,zeta_retro,Theta);

    p_retro_tf = p_collision(x(end),y(end),sigx,sigy,s1+s2);

    x = bplane2x(state_retro_guess(1),state_retro_guess(2),Theta);
    y = bplane2y(state_retro_guess(1),state_retro_guess(2),Theta);
    p_retro(i) = p_collision(x(end),y(end),sigx,sigy,s1+s2);

    % Prograde guess solution
    sol_pro = bvp4c( ...
        @ (t,y) CAM_ode(t,y,kappa,v1,n1,r1,a0,dtheta_t,theta_c), ...
        @ (ya,yb) CAM_bc(ya,yb,r2-r1,0,Q1), ...
        solinit_pro, ...
        options  ...
    );

    t_pro = sol_pro.x*n1/2/pi;
    xi_pro = sol_pro.y(1,:);
    zeta_pro = sol_pro.y(2,:);

    x = bplane2x(xi_pro,zeta_pro,Theta);
    y = bplane2y(xi_pro,zeta_pro,Theta);

    p_pro_tf = p_collision(x(end),y(end),sigx,sigy,s1+s2);

    x = bplane2x(state_pro_guess(1),state_pro_guess(2),Theta);
    y = bplane2y(state_pro_guess(1),state_pro_guess(2),Theta);
    p_pro(i) = p_collision(x(end),y(end),sigx,sigy,s1+s2);

    % Select better of the two guesses
    if p_pro_tf < p_retro_tf
        p_optimal(i) = p_pro_tf;
    else
        p_optimal(i) = p_retro_tf;
    end

end


%% Plots
figure(1);
semilogy(revs,p_optimal,LineWidth=2); grid on; box on; hold on;
semilogy(revs,p_retro,LineWidth=2,LineStyle="--");
semilogy(revs,p_pro,LineWidth=2,LineStyle="--");
xlabel('$n_{revs}$','Interpreter','latex');
ylabel('P_{collision}');


%% Functions
function res = CAM_bc(ya,yb,xi0,zeta0,Q)
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
%         lambda2
%     ];
%---------------------------------------------------------

% final time condition from Pontryagin's maximum principle
lambda_tf = Q*[yb(1);yb(2)];

res = [
    ya(1) - xi0;
    ya(2) - zeta0;
    yb(3) - lambda_tf(1);
    yb(4) - lambda_tf(2);
];

end