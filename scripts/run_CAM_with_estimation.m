%% Notes

% all functions necessary are included at the bottom of this file 
% besides plot_Earth_AAE590

% plot_Earth_AAE590 function and 2k_Earth_July_topography picture
% for generating Earth image (on github)

% satellite measurement info to mess around with:   dt, num
% noise/variance info to mess around with:          sig_r, sig_v, Qs
% orbital parameters to mess around with:           i1, i2, raan1, raan2, h1, h2

% raan1 and raan2 should be equally spaced in order to garuntee TCA at t0

% sig_x and sig_y correspond to rtn frame I believe

% dot(r,u) ~= 0 at all times (tangential CAM is almost optimal)


%% Clear & close

clear;
clc;

tic

format long
fclose('all');
close all

%% Settings
Tolerance = 1e-12;
ode_opts = odeset('RelTol', Tolerance, 'AbsTol', Tolerance);
bvp_opts = bvpset('Stats','off','RelTol',1e-8);


%% Constants
re = 6371e3;  % radius of Earth (m)
mu = 3.986e14; % Earth gravitational parameter (m^3/s^2)


%% Spacecraft physical parameters
s1 = 8;  % Radius of vehicle 1 (m)
s2 = 8;  % Radius of vehicle 2 (m)
sA = s1 + s2;
m1 = 300;
f1 = 3e-3; 
a0 = f1 / m1;  % maximum acceleration magnitude


%% Orbital params at time of closest approach (TCA)

h1 = 550.0e3;  % Orbit height of vehicle 1 (m)
h2 = 550.2e3;  % Orbit height of vehicle 2 (m)

a1 = re + h1;  % Orbital radius of vehicle 1 (m)
a2 = re + h2;  % Orbital radius of vehicle 2 (m)

i1 = 53 * pi/180;  % inclination of S1 (rad)
i2 = 53 * pi/180;  % inclination of S2 (rad)

RAAN1 = 180 * pi/180;  % RAAN of S1 (rad)
RAAN2 =   0 * pi/180;  % RAAN of S2 (rad)

n1 = sqrt(mu/a1^3);  % mean motion of S1 (rad/s)
n2 = sqrt(mu/a2^3);  % mean motion of S2 (rad/s)

T1 = 2*pi/n1;  % period of S1 (s)
T2 = 2*pi/n2;  % period of S2 (s)

theta1 = deg2rad(0); % true anom of S1 at t0 for control maneuver (rad)
theta2 = deg2rad(180); % true anom of S2 at t0 for control maneuver (rad)

% mutual inclination, angle between orbital planes
kappa = acos(sin(i1)*sin(i2)*cos(RAAN1-RAAN2) + cos(i1)*cos(i2));

% Cartesian state at TCA
[r1_vec_true,v1_vec_true] = orb2cartesian(a1,0,i1,RAAN1,0,theta1,mu);
[r2_vec_true,v2_vec_true] = orb2cartesian(a2,0,i2,RAAN2,0,theta2,mu);


%% back-propagate S1 and S2 to get some noisy measurements to use for control problem

t0 = 0;
dt = 20;                                % time step of measurements
num = 50;                              % number of measurements
tk_vec = 0:dt:num*dt-1;                 % measurement time vector

[~,x_S1_backprop] = ode45(@eom_tbp_state, flip(tk_vec), ...
    [r1_vec_true; v1_vec_true], ode_opts, mu);        % backpropagating S1
x_S1_backprop = flip(x_S1_backprop)';

[~,x_S2_backprop] = ode45(@eom_tbp_state, flip(tk_vec), ...
    [r2_vec_true; v2_vec_true], ode_opts, mu);        % backpropagating S2
x_S2_backprop = flip(x_S2_backprop)';


x0_S1_true = x_S1_backprop(:,1);          % saying this is the true state at first measurement time (in actually we wont know this but its nice to show)
x0_S2_true = x_S2_backprop(:,1);
z1 = x_S1_backprop(1:3,:);               % position only measurements
z2 = x_S2_backprop(1:3,:); 

sig_r = 100;                            % position noise (m)
sig_v = 0.1;                            % velocity noise (m/s)
for k = 1:num
    z1(:,k) = z1(1:3,k)+sig_r*randn(3,1); % adding gaussian noise to observations
    z2(:,k) = z2(1:3,k)+sig_r*randn(3,1); % adding gaussian noise to observations
end


%% Perform iterative improvement via recursive leastsquares (LUMVE)

P0bar_S1 = diag(inf(1,6));                 % initial covariance guess is infinity (worst case scenario)
P0bar_S2 = diag(inf(1,6));

x0_S1_ref_start = x0_S1_true+...
    [sig_r*randn(3,1); sig_v*randn(3,1)]; % initial measured state (with noise)
x0_S2_ref_start = x0_S2_true+...
    [sig_r*randn(3,1); sig_v*randn(3,1)];

Rk = diag([sig_r^2,sig_r^2,sig_r^2]);   % measurement noise matrix


% NL RLS minimum variance estimate at first measurement
iter = 3;           % just has to be high enough to converge
show_plots = 1;     % set to 0 for no plots
[x0_S1_ref,P0_S1] = LUMVE(x0_S1_ref_start,P0bar_S1,z1,tk_vec,Rk,mu,iter,show_plots,'Satellite 1',ode_opts);
[x0_S2_ref,P0_S2] = LUMVE(x0_S2_ref_start,P0bar_S1,z2,tk_vec,Rk,mu,iter,show_plots,'Satellite 2',ode_opts);


%% Perform EKF to get covariance and state at TCA

m0_S1 = x0_S1_true;   % initial mean (using for now)
m0_S2 = x0_S2_true;   % initial mean

Qs = (1e-9)^2 * eye(3);    % process noise
M = [zeros(3);             % process noise mapping matrix
    eye(3)];

L = eye(3);     % measurement noise mapping matrix?

% EKF propagation to get state and covariance at TCA
[xkp1,mkp1,Pkp1] = EKF(x0_S1_ref,m0_S1,P0_S1,z1,tk_vec,Rk,Qs,M,L,mu,1,'Satellite 1',ode_opts);
[xkp2,mkp2,Pkp2] = EKF(x0_S2_ref,m0_S2,P0_S2,z2,tk_vec,Rk,Qs,M,L,mu,1,'Satellite 1',ode_opts);

% estimated state of satellite (to use in controller)
r1_vec = xkp1(1:3); 
v1_vec = xkp1(4:6);
r2_vec = xkp2(1:3);
v2_vec = xkp2(4:6);


%% covariance rotation

Cr_S1 = [Pkp1(1,1:3);
         Pkp1(2,1:3);
         Pkp1(3,1:3)];
Cr_S2 = [Pkp2(1,1:3);
         Pkp2(2,1:3);
         Pkp2(3,1:3)];
Cr = Cr_S1 + Cr_S2;

% rotation from ECI to enc frame to find sigmas in enc frame and Theta
[~,DCM_TCA] = eci2enc(r1_vec,v1_vec,r2_vec,v2_vec);
Cr_enc = DCM_TCA*Cr*DCM_TCA';
C = [Cr_enc(1,1) Cr_enc(1,3);
     Cr_enc(3,1) Cr_enc(3,3)];
sig_xi = sqrt(Cr_enc(1,1));
sig_zeta = sqrt(Cr_enc(3,3));
rho_xi_zeta = Cr_enc(3,1)/(sig_xi*sig_zeta);
Theta = 1/2*atan2(2*rho_xi_zeta*sig_xi*sig_zeta,(sig_xi^2-sig_zeta^2));


% rotation from ECI to rtn frame to get sig_x and sig_y 
[Cr_rtn] = eci2rtn(Cr,r1_vec,v1_vec);
sig_x = sqrt(Cr_rtn(1,1));
sig_y = sqrt(Cr_rtn(2,2));
sig_xi_test = sqrt( ...
        (sig_x^2 + sig_y^2 + (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
              );
sig_zeta_test = sqrt( ...
        (sig_x^2 + sig_y^2 - (sig_x^2 - sig_y^2)*cos(2*Theta)) / 2 ...
                );
rho_xi_zeta_test = (sig_x^2 - sig_y^2) * sin(2*Theta) / sig_xi / sig_zeta;


Q = [
        1/sig_xi^2    -rho_xi_zeta/sig_xi/sig_zeta;
        -rho_xi_zeta/sig_xi/sig_zeta  1/sig_zeta^2;
    ];


%% Timing paramters (thrust arc and coast arc)

npts = 10000;
del_theta_t = 180 * pi/180; % thrust arc duration (mess around w this)
n_revs = pi - del_theta_t/2; % next conjunction 180 degrees after first
del_theta_c = 2*pi*n_revs - del_theta_t/2;
theta_c = 2*pi*n_revs + del_theta_t/2;

tf = del_theta_t/n1+t0;
t = linspace(t0,tf,npts);


%% TPBVP Initial conditions and necessary parameters

r1 = norm(r1_vec);
r2 = norm(r2_vec);
v1 = norm(v1_vec);
v2 = norm(v2_vec);

[r1_vec_rtn] = eci2rtn(r1_vec,r1_vec,v1_vec);
[v1_vec_rtn] = eci2rtn(v1_vec,r1_vec,v1_vec);
[r2_vec_rtn] = eci2rtn(r2_vec,r1_vec,v1_vec);
[v2_vec_rtn] = eci2rtn(v2_vec,r1_vec,v1_vec);
[r1_vec_enc] = eci2enc(r1_vec,v1_vec,r2_vec,v2_vec);

xi0 = r1_vec_enc(1);
zeta0 = r1_vec_enc(3);
b0 = [xi0; zeta0];


%% Run TPBVP

% obtaining lambda0 guess using Q*bf using tangential CAM eqns
del_r1 = 2*a0*r1^3/mu*(del_theta_t+sin(del_theta_c)-sin(del_theta_t+del_theta_c));
del_t = a0*r1^(7/2)/mu^(3/2)*(3*del_theta_t*(del_theta_t/2+del_theta_c)- ...
        8*sin(del_theta_t/2)*sin(del_theta_t/2+del_theta_c));
if r2 > r1
    xif_guess = xi0 - del_r1;
elseif r1 > r2
    xif_guess = xi0 + del_r1;
end
zetaf_guess = zeta0+a0*r1^3/mu*cos(kappa/2)* ...
                (3*del_theta_t*(del_theta_t/2+del_theta_c)- ...
                8*sin(del_theta_t/2)*sin(del_theta_t/2+del_theta_c));
bf_guess = [xif_guess; zetaf_guess];
lambda0_guess = Q*bf_guess;


solinit = bvpinit(t, [b0; lambda0_guess]);
sol = bvp4c( ...
    @ (t,y) bvp_ode(t,y,kappa,v1,n1,r1,a0,theta1,theta_c,r2), ...
    @ (ya,yb) bvp_bc(ya,yb,b0,Q), ...
    solinit, ...
    bvp_opts  ...
);

xi = sol.y(1,:);
zeta = sol.y(2,:);
b = sol.y(1:2,:);

lambda1 = sol.y(3,:);
lambda2 = sol.y(4,:);
lambda = sol.y(3:4,:);


bplane2x = @ (xi,zeta,Theta) xi*cos(Theta) + zeta*sin(Theta);
bplane2y = @ (xi,zeta,Theta) xi*sin(Theta) - zeta*cos(Theta);

x = bplane2x(xi,zeta,Theta);
y = bplane2y(xi,zeta,Theta);

p_t0_Ayuso = p_collision_calc_Ayuso(x(1),y(1),sig_x,sig_y,sA)
p_tf_Ayuso = p_collision_calc_Ayuso(x(end),y(end),sig_x,sig_y,sA)

% trunc = 4;
% p_t0_Chan = p_collision_calc_Chan(xi(1),zeta(1),sig_xi,sig_zeta,rho_xi_zeta,sA,trunc);
% p_tf_Chan = p_collision_calc_Chan(xi(end),zeta(end),sig_xi,sig_zeta,rho_xi_zeta,sA,trunc);


%% finding control input and propagating S1 with control input

u = NaN(3,length(t));
tangential_test = NaN(1,length(t));
dbdt = NaN(2,length(t));
%p_Chan = NaN(1,length(t));
p_Ayuso = NaN(1,length(t));

% control input calc
for ct=1:length(t)

    theta = theta1 + t(ct) * n1;
    dtheta = theta_c - theta;
    
    R = [
                  0.0,          0.0,       -1.0;    
            -cos(kappa/2), -sin(kappa/2),   0.0
        ];
    
    if r2 > r1
        K = [
                -v1,         0.0,       0.0;
                 0.0,        0.0,       1.0;
                 0.0,        1.0,       0.0
            ];
    elseif r2 < r1
        K = [
                -v1,         0.0,       0.0;
                 0.0,        0.0,      -1.0;
                 0.0,       -1.0,       0.0
            ]; 
    end 
    
    Drr  = sin(dtheta) / n1;
    Drth = (2 - 2*cos(dtheta)) / n1;
    Dtr  = (2 - 2*cos(dtheta)) / (n1^2 * a1);
    Dtth = (3*dtheta - sin(dtheta)) / (n1^2 * a1);  % Works without the 4 in front?
    Dwh  = sin(dtheta) / n1;
    
    D = [
            Dtr, Dtth, 0.0;
            Drr, Drth, 0.0;
            0.0,  0.0, Dwh
        ];
    
    M = R*K*D;
    
    u(:,ct) = a0*M'*lambda(:,ct) / norm(M'*lambda(:,ct));
    dbdt(:,ct) = M * u(:,ct);
    p_Ayuso(ct) = p_collision_calc_Ayuso(x(ct),y(ct),sig_x,sig_y,sA);
    %p_Chan(ct) = p_collision_calc_Chan(xi(ct),zeta(ct),sig_xi,sig_zeta,rho_xi_zeta,sA,trunc);

end
u_hat = u./vecnorm(u);


% propagating S1 w/ control input and uncontrolled S2
% may want to plot state vs time cus its hard to see with a small u_max
[~,x_S1_controlled] = ode45(@eom_tbp_state_controlled, t ,[r1_vec; v1_vec], ode_opts, mu, t, u);  
x_S1_controlled = x_S1_controlled';
[~,x_S2] = ode45(@eom_tbp_state, t ,[r2_vec; v2_vec], ode_opts, mu);  
x_S2 = x_S2';

for k=1:npts
    r1_rtn = eci2rtn(x_S1_controlled(1:3,k),x_S1_controlled(1:3,k),x_S1_controlled(4:6,k));
    tangential_test(k) = dot(r1_rtn,u(:,k))/(norm(r1_rtn)*norm(u(:,k)));
end


% plot of measurements/path before TCA and controlled path after TCA
figure()
scatter3(z1(1,:)/1000,z1(2,:)/1000,z1(3,:)/1000,5,'filled','r')
hold on
plot3(x_S1_backprop(1,:)/1000,x_S1_backprop(2,:)/1000,x_S1_backprop(3,:)/1000,'-r')
hold on
scatter3(z2(1,:)/1000,z2(2,:)/1000,z2(3,:)/1000,5,'filled','g')
hold on
plot3(x_S2_backprop(1,:)/1000,x_S2_backprop(2,:)/1000,x_S2_backprop(3,:)/1000,'-g')
hold on
plot3(x_S1_controlled(1,:)/1000,x_S1_controlled(2,:)/1000,x_S1_controlled(3,:)/1000,'--r')
hold on
plot3(x_S2(1,:)/1000,x_S2(2,:)/1000,x_S2(3,:)/1000,'-g')
hold on
plot_Earth_AAE590;
grid on
title('$S_1$ and $S_2$ pre/post TCA','Interpreter','latex')
legend('$S_1$ measurements pre TCA','$S_1$ true path pre TCA', ...
    '$S_2$ measurements pre TCA','$S_2$ true path pre TCA', ...
    'controlled $S_1$ path post TCA', '$S_2$ path post TCA', ...
    'interpreter','latex')


% plot of state and control history
figure()
subplot(2,2,1);
plot(t,b(1,:))
hold on
plot(t,b(2,:))
grid on
title('$b$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$b \; [m]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$\xi$','$\zeta$','Interpreter','latex')
subplot(2,2,2);
plot(t,dbdt(1,:))
hold on
plot(t,dbdt(2,:))
grid on
title('$\dot{b}$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$\dot{b} \; [\frac{m}{s}]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$\dot{\xi}$','$\dot{\zeta}$','Interpreter','latex')
subplot(2,2,3);
plot(t,u_hat(1,:))
hold on
plot(t,u_hat(2,:))
hold on
plot(t,u_hat(3,:))
grid on
title('$u$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$u \; [\frac{m}{s^2}]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$u_{r}$','$u_{\theta}$','$u_h$','Interpreter','latex')
subplot(2,2,4);
plot(t,tangential_test)
grid on
title('$\hat{u} \cdot \hat{r}$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$\hat{u} \cdot \hat{r}$','Interpreter','latex')
xlim([t(1) t(end)])


% plot of state history in ECI frame
figure()
subplot(2,2,1);
plot(t,x_S1_controlled(1,:),'r')
hold on
plot(t,x_S2(1,:),'g')
grid on
title('$x$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$x \; [m]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$x_{S_1}$','$x_{S_2}$','Interpreter','latex')
subplot(2,2,2);
plot(t,x_S1_controlled(2,:),'r')
hold on
plot(t,x_S2(2,:),'g')
grid on
title('$y$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$y \; [m]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$y_{S_1}$','$y_{S_2}$','Interpreter','latex')
subplot(2,2,3);
plot(t,x_S1_controlled(1,:),'r')
hold on
plot(t,x_S2(1,:),'g')
grid on
title('$z$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$z \; [m]$','Interpreter','latex')
xlim([t(1) t(end)])
legend('$z_{S_1}$','$z_{S_2}$','Interpreter','latex')
subplot(2,2,4);
plot(t,p_Ayuso,'k')
grid on
title('$P_{collision}$ vs $t$','Interpreter','latex')
xlabel('$t \; [s]$','Interpreter','latex')
ylabel('$P_{collision}$','Interpreter','latex')
xlim([t(1) t(end)])


toc


%% Functions  

% TPBVP functions
function res = bvp_bc(ya,yb,b0,Q)

    % final time condition from Pontryagin's maximum principle
    lambda_tf = Q*yb(3:4);

    res =   [
                ya(1:2) - b0;
                yb(3:4) - lambda_tf;
            ];

end

function dydt = bvp_ode(t,y,k,v1,n1,r1,a0,theta10,theta_c,r2)
    
    b = y(1:2);
    lambda = y(3:4);
    
    % Parameters for state update and input updates
    theta = theta10 + t * n1;
    dtheta = theta_c - theta;
    
    % Optimal Input from Pontryagin Maximum Principle
    R = [
                  0.0,       0.0,  -1.0;
            -cos(k/2), -sin(k/2),   0.0
        ];
    
    if r2 > r1
        K = [
                -v1,         0.0,       0.0;
                 0.0,        0.0,       1.0;
                 0.0,        1.0,       0.0
            ];
    elseif r2 < r1
        K = [
                -v1,         0.0,       0.0;
                 0.0,        0.0,      -1.0;
                 0.0,       -1.0,       0.0
            ]; 
    end

    Drr  = sin(dtheta) / n1;
    Drth = (2 - 2*cos(dtheta)) / n1;
    Dtr  = (2 - 2*cos(dtheta)) / (n1^2 * r1);
    Dtth = (3*dtheta - sin(dtheta)) / (n1^2 * r1);  % Works without the 4 in front?
    Dwh  = sin(dtheta) / n1;
    
    D = [
            Dtr, Dtth, 0.0;
            Drr, Drth, 0.0;
            0.0,  0.0, Dwh
        ];
    
    M = R*K*D;
    
    u = a0*M'*lambda / norm(M'*lambda);
    
    dbdt = M * u;
    dlambdadt = 0;
    
    dydt(1:2) = dbdt;
    dydt(3:4) = dlambdadt;

end


% reference frame rotation functions
function [array_rtn,DCM] = eci2rtn(array_rot,r_vec_ref,v_vec_ref)

    % correcting dimensions
    r_vec_ref = reshape(r_vec_ref,3,1);
    v_vec_ref = reshape(v_vec_ref,3,1);

    % defining unit vectors
    r_hat = r_vec_ref./norm(r_vec_ref);
    n_hat = cross(r_vec_ref,v_vec_ref)./norm(cross(r_vec_ref,v_vec_ref));
    t_hat = cross(n_hat,r_hat);

    % rotation
    DCM = [r_hat t_hat n_hat]';

    if length(array_rot(:)) == 3 % vector rotation
        array_rtn = DCM*array_rot;
    elseif length(array_rot(:)) == 9 % matrix rotation
        array_rtn = DCM*array_rot*DCM';
    end

end

function [r1_vec_enc,DCM] = eci2enc(r1_vec,v1_vec,r2_vec,v2_vec)

    % correcting dimensions
    r1_vec = reshape(r1_vec,3,1);
    v1_vec = reshape(v1_vec,3,1);
    r2_vec = reshape(r2_vec,3,1);
    v2_vec = reshape(v2_vec,3,1);

    % del_r centered at S2
    del_r_vec = r1_vec-r2_vec;

    % defining unit vectors
    xi_vec = cross(v2_vec,v1_vec);
    xi_hat = xi_vec./norm(xi_vec);
    eta_vec = v1_vec-v2_vec;
    eta_hat = eta_vec./norm(eta_vec);
    zeta_hat = cross(xi_hat,eta_hat);
    
    % rotation
    DCM = [xi_hat eta_hat zeta_hat]';
    r1_vec_enc = DCM*del_r_vec;  

end


% estimation functions
function [x0_ref,P0] = LUMVE(x0_ref_start,P0bar,z,t_vec,Rk,mu,iter,flag,satname,ode_opts)

    x0_ref = x0_ref_start;
    t0 = t_vec(1);
    
    delx0bar = zeros(6,1);

    Ri = Rk^-1;                             % inverse of Rk
    
    % begin the iteration loop
    for loop = 1:iter
        
        % initialize time variable, the reference state, and the STM
        Lam = inv(P0bar);
        lam = P0bar\delx0bar;
        tkm1 = t0;
        xref = x0_ref;
        STM = eye(6);
    
        % initialize storage for the residual and the time of the residual
        rest = zeros(1,length(t_vec));
        resm = zeros(3,length(t_vec));
    
        for k = 1:length(t_vec)
    
            % extract the time, measurement, covariance
            tk = t_vec(k);
            zk = z(:,k);
    
            % propagate the ref state and STM, but only if were past t0
            if tk > t0
                [~,XX] = ode45(@eom_tbp_state, [tkm1, tk], [xref; STM(:)], ode_opts, mu);
                xref = XX(end, 1:6)';
                STM = reshape(XX(end, 7:end)', 6, 6);
            end
    
            % compute ref measurement
            zref = xref(1:3);
    
            % measurement mapping matrix (Htilde)
            Htilde = [eye(3) zeros(3)];
            
            % accumulate the lambdas for LUMVE
            H = Htilde*STM;
            lam = lam + H'*Ri*(zk - zref); 
            Lam = Lam + H'*Ri*H;
            
            % store the time and the residual for plotting later
            resm(:,k) = zk - zref;
            rest(k) = tk;
    
            % reset the timing variable
            tkm1 = tk;
        end
    
        % get the least squares solution
        delx0hat = Lam\lam;
        P0 = Lam^-1;
        
        % perform the iteration by shifting the reference and the estimated deviation
        x0_ref = x0_ref + delx0hat;
        delx0bar = delx0bar - delx0hat;
    
        % plotting residuals if wanted (shows convergence)
        if flag 
            titleString = sprintf('%s: $r$ residual over time for iteration %d', satname, loop);
            figure()
            scatter(rest,resm(1,:),15,'filled','r')
            hold on 
            scatter(rest,resm(2,:),15,'filled','b')
            hold on 
            scatter(rest,resm(3,:),15,'filled','g')
            grid on
            title(titleString,'Interpreter','latex');
            xlabel('$t \; [s]$','Interpreter','latex')
            ylabel('$\epsilon_r \; [m]$','Interpreter','latex')
            legend('$\epsilon_x$','$\epsilon_y$','$\epsilon_z$','interpreter','latex')
        end
    
    end
end

function [xk,mkp,Pkp] = EKF(x0,m0,P0,z,t_vec,Rk,Qs,M,L,mu,flag,satname,ode_opts)
    
    nx = length(x0); % state size
    [nq,~] = size(Qs);
    [nr,~] = size(L);

    t0 = t_vec(1);
    
    % storage for err and stdev for plotting
    tplot = zeros(1,2*length(t_vec)-1);            % time (s)
    eplot = zeros(nx, 2*length(t_vec)-1);          % state error (m; m/s)
    splot = zeros(nx, 2*length(t_vec)-1);          % state stdev (m; m/s)
    rplot = zeros(nr, length(t_vec)-1);            % measurement residuals (m)
    wplot = zeros(nr, length(t_vec)-1);            % measurement stdev (m)
    ctr = 1;                                        % counter variable for loop
    tplot(:, ctr) = t0;          
    eplot(:, ctr) = x0 - m0;          
    splot(:, ctr) = sqrt(diag(P0));          

    tkm1 = t0;      % initial time (iterated in loop)
    xkm1 = x0;      % initial state (iterated in loop)
    mkm1 = m0;      % initial mean (iterated in loop)
    Pkm1 = P0;      % initial covariance (iterated in loop)

    
    % begin the time loop
    for k = 2:length(t_vec)
    
        tk = t_vec(k);
    
        % propagate true state (with noise) from tkm1 to tk
        wkm1 = chol(Qs)' * randn(nq, 1);
        [~, XX] = ode45(@eom_tbp_state_noisy,[tkm1, tk],xkm1,ode_opts,mu,M,wkm1);
        xk = XX(end, :)';
    
        % propagate mean and covariance from time tkm1 to tk
        [~,XX] = ode45(@eom_tbp_meancov, [tkm1, tk], [mkm1; Pkm1(:)], ode_opts, mu, Qs, M);
        mkm = XX(end, 1:nx)';
        Pkm = reshape(XX(end, nx + 1:end)', nx, nx);
        Pkm = 0.5 * (Pkm + Pkm'); % assures Pkm is symmetric
        
        % measurement and ref measurement
        zk = z(:,k);
        zhatk = xk(1:3);
    
        % update mean and covariance at time tk
        Hk = [eye(3) zeros(3)];
        Ck = Pkm*Hk';
        Wk = Hk*Pkm*Hk' + L*Rk*L';
        Kk = Ck/Wk;
        mkp = mkm + Kk*(zk - zhatk);
        Pkp = Pkm - Kk*Hk*Pkm;
        Pkp = 0.5 * (Pkp + Pkp');
    
        % store the a prior / a posteriori error and stdev
        ctr = ctr + 1;
        tplot(:, ctr) = tk;
        eplot(:, ctr) = xk - mkm;
        splot(:, ctr) = sqrt(diag(Pkm));
        ctr = ctr + 1;
        tplot(:, ctr) = tk;
        eplot(:, ctr) = xk - mkp;
        splot(:, ctr) = sqrt(diag(Pkp));
    
        % store the residual and its standard deviation
        rplot(:, k-1) = zk - zhatk;
        wplot(1, k-1) = sqrt(Wk(1,1));
        wplot(2, k-1) = sqrt(Wk(2,2));
        wplot(3, k-1) = sqrt(Wk(3,3));
    
        % cycle the time, state, mean, and covariance
        tkm1 = tk;
        xkm1 = xk;
        mkm1 = mkp;
        Pkm1 = Pkp;
    
    end
    
    if flag
        % plot the estimation errors and 3-sigma "bounds"
        xlab = {'$t \; [s]$'};
        ylab = {'$\delta_x \; [m]$', '$\delta_y \; [m]$', '$\delta_z \; [m]$', '$\delta_{v_x} \; [\frac{m}{s}]$', '$\delta_{v_y} \; [\frac{m}{s}]$', '$\delta_{v_z} \; [\frac{m}{s}]$'};
        tit = {'$x$-position error vs time', '$y$-position error vs time','$z$-position error vs time','$x$-velocity error vs time','$y$-velocity error vs time','$z$-velocity error vs time'};
        figure()
        for i =1:nx
            subplot(2,3,i);
            C = get(gca, 'ColorOrder');
            plot(tplot(1,:), eplot(i,:), 'Color', C(1,:), 'LineWidth', 1.2);
            hold on
            plot(tplot(1,:), 3*splot(i,:), 'Color', C(2,:), 'LineWidth', 1.2);
            plot(tplot(1,:), -3*splot(i,:), 'Color', C(2,:), 'LineWidth', 1.2);
            xlabel(xlab,'interpreter','latex');
            ylabel(ylab{i},'interpreter','latex');
            title(tit(i),'interpreter','latex');
            legend('errors','3-sigma bound','interpreter','latex')
        end
        titleString = sprintf('For %s:', satname);
        sgtitle(titleString,'interpreter','latex')
        
        % plot the measurement residuals and 3-sigma "bounds"
        figure()
        subplot(3,1,1);
        C = get(gca, 'ColorOrder');
        plot(t_vec(2:end), rplot(1,:), 'x', 'Color', C(1,:))
        hold on
        plot(t_vec(2:end), 3*wplot(1,:), 'Color', C(2,:), 'LineWidth', 1.2);
        hold on
        plot(t_vec(2:end), -3*wplot(1,:), 'Color', C(2,:), 'LineWidth', 1.2);
        title('$x$ measurement residuals vs time','interpreter','latex')
        xlabel('$t \; [s]$','interpreter','latex')
        ylabel('$\epsilon_x \; [m]$','interpreter','latex')
        xlim([t_vec(1) t_vec(end)+5])
        legend('residuals','3-sigma bound','interpreter','latex')
        subplot(3,1,2);
        C = get(gca, 'ColorOrder');
        plot(t_vec(2:end), rplot(2,:), 'x', 'Color', C(1,:))
        hold on
        plot(t_vec(2:end), 3*wplot(2,:), 'Color', C(2,:), 'LineWidth', 1.2);
        hold on
        plot(t_vec(2:end), -3*wplot(2,:), 'Color', C(2,:), 'LineWidth', 1.2);
        title('$y$ measurement residuals','interpreter','latex')
        xlabel('$t \; [s]$','interpreter','latex')
        ylabel('$\epsilon_y \; [m]$','interpreter','latex')
        xlim([t_vec(1) t_vec(end)+5])
        legend('residuals','3-sigma bound','interpreter','latex')
        subplot(3,1,3);
        C = get(gca, 'ColorOrder');
        plot(t_vec(2:end), rplot(3,:), 'x', 'Color', C(1,:))
        hold on
        plot(t_vec(2:end), 3*wplot(3,:), 'Color', C(2,:), 'LineWidth', 1.2);
        hold on
        plot(t_vec(2:end), -3*wplot(3,:), 'Color', C(2,:), 'LineWidth', 1.2);
        title('$z$ measurement residuals','interpreter','latex')
        xlabel('$t \; [s]$','interpreter','latex')
        ylabel('$\epsilon_z \; [m]$','interpreter','latex')
        xlim([t_vec(1) t_vec(end)+5])
        legend('residuals','3-sigma bound','interpreter','latex')

        sgtitle(titleString,'interpreter','latex')
    end
end


% two body problem equations of motion functions
function [dx] = eom_tbp_state_controlled(t, x, mu, t_hist, u_hist)
    
    % Identify the position and velocity from the state vector:
    r_vec = x(1:3); %[km]
    v_vec = x(4:6); %[km s-1]
    r = norm(r_vec);

    % control interpolation
    u1_rtn = interp1(t_hist,u_hist(1,:),t);
    u2_rtn = interp1(t_hist,u_hist(2,:),t);
    u3_rtn = interp1(t_hist,u_hist(3,:),t);
    u_rtn = [u1_rtn; u2_rtn; u3_rtn];

    % control rotation
    r_hat = r_vec/r;
    n_hat = cross(r_vec,v_vec)/norm(cross(r_vec,v_vec));
    t_hat = cross(n_hat,r_hat);
    DCM = [r_hat t_hat n_hat]';
    u_eci = DCM^-1*u_rtn;

    % Assign values in the state's derivative:
    dx = [v_vec; -mu*r_vec./r^3 + u_eci];

end

function[dx] = eom_tbp_state(t,x,mu)

    % Identify the position and velocity from the state vector:
    r = x(1:3); %[km]
    v = x(4:6); %[km s-1]
    
    % Assign values in the state's derivative:
    dx = zeros(6,1);
    dx(1:3) = v; %[km s-1] velocity
    dx(4:6) = -mu/norm(r)^3.*r; %[km s-2] acceleration
    
    if length(x) > 6 
        Phi = reshape(x(7:end), 6, 6);
        Gr = mu/norm(r)^5*(3*(r*r')-norm(r)^2*eye(3));
        F = [zeros(3) eye(3);
             Gr zeros(3)];
        dPhi = F*Phi;
        dx = [dx; dPhi(:)];
    end
end

function[dx] = eom_tbp_state_noisy(t,x,mu,M,w)

    % Identify the position and velocity from the state vector:
    r = x(1:3); %[km]
    v = x(4:6); %[km s-1]
    
    % Assign values in the state's derivative:
    dx = zeros(6,1);
    dx(1:3) = v; %[km s-1] velocity
    dx(4:6) = -mu/norm(r)^3.*r; %[km s-2] acceleration
    dx(1:6) = dx(1:6)+M*w;

end

function[dx] = eom_tbp_meancov(t,x,mu,Qs,M)

    % Identify the position and velocity from the state vector:
    r = x(1:3); %[km]
    v = x(4:6); %[km s-1]
    
    % Assign values in the state's derivative:
    dx = zeros(6,1);
    dx(1:3) = v; %[km s-1] velocity
    dx(4:6) = -mu/norm(r)^3.*r; %[km s-2] acceleration
    
    if length(x) > 6
        P = reshape(x(7:end), 6, 6);
        Gr = mu/norm(r)^5*(3*(r*r')-norm(r)^2*eye(3));
        F = [zeros(3) eye(3);
             Gr zeros(3)];
        dP = F*P + P*F' + M*Qs*M';
        dx = [dx; dP(:)];
    end
end


% orbital elements to Cartesian state functions
function [r_vec,v_vec] = orb2cartesian(a,e,inc,raan,argp,f,mu)
    
    p = a*(1-e^2); % semi-latus rectum [k]

    r_pi = [cos(f); sin(f); 0];
    v_pi = [-sin(f); e+cos(f); 0];

    DCM = [cos(raan)*cos(argp)-sin(raan)*sin(argp)*cos(inc), -cos(raan)*sin(argp)-sin(raan)*cos(argp)*cos(inc), sin(raan)*sin(inc); ...
        sin(raan)*cos(argp)+cos(raan)*sin(argp)*cos(inc), -sin(raan)*sin(argp)+cos(raan)*cos(argp)*cos(inc), -cos(raan)*sin(inc); ...
        sin(argp)*sin(inc), cos(argp)*sin(inc), cos(inc)];

    r_vec = DCM*r_pi;
    v_vec = DCM*v_pi;
    r_vec = p*r_vec/(1+e*cos(f));
    v_vec = sqrt(mu/p)*v_vec;

end


% collision probability functions
function P = p_collision_calc_Ayuso(x,y,sig_x,sig_y,S_A)

    P = S_A^2 / (2*sig_x*sig_y) * ...
        (1 + 1/8*((x.^2/sig_x^2 - 1)/(sig_x^2) + ((y.^2/sig_y^2 - 1)/(sig_y^2) )) ) .* ...
        exp(-1/2*(x.^2/sig_x^2 + y.^2/sig_y^2));

end

% smth is wrong with this function below
function P = p_collision_calc_Chan(xi,zeta,sig_xi,sig_zeta,rho_xi_zeta,sA,trunc)
    
    u = sA^2/(sig_xi*sig_zeta*sqrt(1-rho_xi_zeta^2));
    v = ((xi/sig_xi)^2+(zeta/sig_zeta)^2-2*rho_xi_zeta*xi/sig_xi*zeta/sig_zeta)/(1-rho_xi_zeta^2);

    ksum = 0;
    msum = 0;
    for m = 0:trunc
        for k = 0:m
            ksum = ksum + u^k/(2^k*factorial(k));
        end
        msum = msum + v^m/(2^m*factorial(m))*(1-exp(-u/2)*ksum);
    end
    P = exp(-v/2)*msum;

end
