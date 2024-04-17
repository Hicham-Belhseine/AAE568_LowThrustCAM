function dydt = CAM_ode(y,t,k,v,u,a0,Q,dtheta_t)
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
%     k - angle between the two orbital planes (mutual inclination)
%     v - Orbital velocity of S1
%     u - argument of latitude
%
% Authors : 
%    Matthew Stephens, steph185@purdue.edu
%    Ryan Herberg, herberg@purdue.edu
%    Hicham Belhseine, hbelhsei@purdue.edu
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

% 


end