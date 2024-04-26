function [costate_retro, costate_pro, state_retro, state_pro] = guess_CAM_costate(n_revs, dtheta_t, a0, r1, kappa, Q, xi_0, zeta_0)
%% Function Name: guess_CAM_costate
%
% Description: Provides a costate guess based on the retrograde tangential
% burn and the prograde tangential burn. 
%
% Assumptions: Circular orbit, negligible mass loss during burn
%
% Inputs:
%     n_revs - number of revs from center of thrust to conjunction event
%     dtheta_t - thrust arc length (radians)
%     a0 - max acceleration of maneuverable satellite
%     r1 - radius of maneuverable satellite orbit
%     kappa - angle between the two orbital planes (mutual inclination)
%     Q - Q matrix in objective function
%     xi_0 - initial condition for xi
%     zeta_0 - initial condition for zeta
%     
%
% Authors : 
%    Matthew Stephens, steph185@purdue.edu
%    Ryan Herberg, herberg@purdue.edu
%    Hicham Belhseine, hbelhsei@purdue.edu
% 
% Date: April 24, 2024
%---------------------------------------------------------
mu = 3.986e14; % Earth gravitational parameter (m^3/s^2)

dtheta_c = 2*pi*n_revs - dtheta_t/2;

xi_retro_tf = xi_0 - 2*a0*(r1^3)/mu * ...
    (dtheta_t + sin(dtheta_c) - sin(dtheta_t+dtheta_c));

zeta_retro_tf = zeta_0 + a0*(r1^3)/mu*cos(kappa/2)* ...
    (3*dtheta_t*(dtheta_t/2 + dtheta_c) - 8*sin(dtheta_t/2)*sin(dtheta_t/2 + dtheta_c));

xi_pro_tf = xi_0 + 2*a0*(r1^3)/mu * ...
    (dtheta_t + sin(dtheta_c) - sin(dtheta_t+dtheta_c));

zeta_pro_tf = zeta_0 - a0*(r1^3)/mu*cos(kappa/2)* ...
    (3*dtheta_t*(dtheta_t/2 + dtheta_c) - 8*sin(dtheta_t/2)*sin(dtheta_t/2 + dtheta_c));

state_retro = [xi_retro_tf; zeta_retro_tf];

state_pro = [xi_pro_tf; zeta_pro_tf];

costate_retro = Q * [xi_retro_tf; zeta_retro_tf];

costate_pro = Q * [xi_pro_tf; zeta_pro_tf];

end