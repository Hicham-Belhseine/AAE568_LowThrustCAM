

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