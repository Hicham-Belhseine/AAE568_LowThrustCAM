function dydt = CAM_ode(y,t,tc,k,v1,n1,a0,dtheta_t)
%% Function Name: CAM_ode
%
% Description: Set of ODEs describing a low-thrust collision avoidance
% maneuver (CAM)
%
% Assumptions: Circular orbit, negligible mass loss
%
% Inputs:
%     y - state vector
%     t - time
%     tc - time of collision
%     k - angle between the two orbital planes (mutual inclination)
%     v1 - Orbital velocity of S1
%     n1 - mean motion
%     
%
% Authors : 
%    Matthew Stephens, steph185@purdue.edu
%    Ryan Herberg, herberg@purdue.edu
%    Hicham Belhseine, hbelhsei@purdue.edu
% 
% Date: April 16, 2024
%---------------------------------------------------------

% 2D position vector of S1 relative to S2 in B-plane
% xi: cross-product of S1/S2 velocity
% eta: along velocity of S1 relative to S2, always zero
% zeta: completes right-handed coordinate system
xi   = y(1); 
zeta = y(2); 

% Lagrangian multiplier/co-state variables
% Constant over time
lambda1 = y(3);
lambda2 = y(4);

% Parameters for state update and input updates
dtheta = (tc - t) / n1;
theta = t / n1;

% Optimal Input from Pontryagin Maximum Principle
R = [
          0.0,       0.0,  -1.0;
    -cos(k/2), -sin(k/2),   0.0
];

K = [
    -v1,         0.0,       0.0;
     0.0,        0.0, -sign(xi);
     0.0,  -sign(xi),       0.0
];

Drr  = sin(dtheta) / n1;
Drth = (2 - 2*cos(dtheta)) / n1;
Dtr  = (2 - 2*cos(dtheta)) / (n1^2 * r1);
Dtth = 4 * (3*dtheta - sin(dtheta)) / (n1^2 * r1);
Dwh  = sin(dtheta) / n1;

D = [
    Dtr, Dtth, 0.0;
    Drr, Drth, 0.0;
    0.0,  0.0, Dwh
];

M = R*K*D;

u = a0*M'*[lambda1; lambda2] / norm(a0*M'*[lambda1; lambda2]);

% apply state update if still part of thrust arc
if theta < dtheta_t
    dbdt = M * u;
else
    dbdt = [0; 0];
end

dydt(1) = dbdt(1);
dydt(2) = dbdt(2);
dydt(3) = 0;
dydt(4) = 0;

end
