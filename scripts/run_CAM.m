%% Scenario Description
% Section VII Scenario in Hernando-Ayuso/Bombardelli paper

%% Parameters
% Altitude/velocity
re = 6371e3;

h1 = 550.0e3;
h2 = 550.2e3; 

r1 = re + h1; 
r2 = re + h2;  

mu = 3.986e14;
v1 = sqrt(mu / r1^3);

% Orbital params
i1 = 53 * pi/180;  % inclination of S1
i2 = 53 * pi/180;  % inclination of S2

RAAN1 = 180 * pi/180;  % RAAN of S1
RAAN2 =   0 * pi/180;  % RAAN of S2

% mutual inclination, angle between orbital planes
kappa = acos(sin(i1)*sin(i2)*cos(180) + cos(i1)*cos(i2));  

% Timing
T1 = 2*pi*sqrt(r1^3 / mu);
T2 = 2*pi*sqrt(r2^3 / mu);

dtheta_t = 200 * pi;

solinit = bvpinit(linspace(0,1,50), [0 0 1 1 1]);
options = bvpset('Stats','on','RelTol',1e-1);

sol = bvp4c(@CAM_ode, @CAM_bc, solinit, options);

%% Functions
function res = CAM_bc(ya,yb)
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

    res = [

    ];

end