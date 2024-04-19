function P = p_collision(x,y,sig_x,sig_y,S_A)
%% Function Name: p_collision
%
% Description: Equation [5]
% 
% Inputs:
%     x - 
%     y - 
%     sig_x - 
%     sig_y - 
%     S_A - 
%     
% Authors : 
%    Matthew Stephens, steph185@purdue.edu
%    Ryan Herberg, herberg@purdue.edu
%    Hicham Belhseine, hbelhsei@purdue.edu
% 
% Date: April 18, 2024
%---------------------------------------------------------

P = S_A^2 / (2*sig_x*sig_y) * ...
    (1 + 1/8*((x^2/sig_x^2 - 1)/(sig_x^2) + ((y^2/sig_y^2 - 1)/(sig_y^2) )) ) * ...
    exp(-1/2*(x^2/sig_x^2 + y^2/sig_y^2));

end
